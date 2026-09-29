/*
 * Copyright (C) 2014-2021 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * Beagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beagle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package phase;

import blbutil.FloatArray;
import blbutil.FloatList;
import ints.IntArray;
import ints.IntList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;
import vcf.Markers;

/**
 * <p>Class {@code PhaseBaum1} implements the forward and backward algorithms
 * for a haploid Li and Stephens hidden Markov model.
 * </p>
 * <p>Instances of class {@code PhaseBaum1} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PhaseBaum1 {

    private final PhaseData phaseData;
    private final boolean burnin;
    private final EstPhase estPhase;
    private final Markers markers;
    private final int nMarkers;
    private final List<int[]> refAlleles;
    private final byte[][][] mismatches;  // [3][nMarker][nStates]
    private final float pMismatch;
    private final float[] emProbs;

    private final int maxStates;
    private final BasicPhaseStates states;
    private final FloatList lrList;

    private int nStates;
    private final float[][] fwd;
    private final float[][] bwd;
    private final float[] sum;
    private final List<float[]> missProbs1;
    private final List<float[]> missProbs2;
    private final List<float[]> bwdHet1;
    private final List<float[]> bwdHet2;

    private boolean swapHaps = false;
    private int missIndex = -1;

    private int swapCnt = 0;
    private static final AtomicLong nSwaps = new AtomicLong(0);
    private static final AtomicLong nUnphHets = new AtomicLong(0);

    /**
     * Returns the proportion of unphased heterozygotes whose phase
     * relative to the previous heterozygote has been changed.
     * The counters for the number of heterozygotes whose phase has been
     * changed and for the total number of heterozygotes are then
     * re-initialized to 0.
     * @return the proportion of unphased heterozygotes whose phase
     * has been changed
     */
    public static double getAndResetSwapRate() {
        double rate = (double) nSwaps.get() / nUnphHets.get();
        nSwaps.set(0);
        nUnphHets.set(0);
        return rate;
    }

    /**
     * Creates a {@code PhaseLSBaum} instance from the specified data.
     *
     * @param phaseIbs the IBS haplotype segments

     * @throws NullPointerException if {@code phaseIBS == null}
     */
    public PhaseBaum1(PbwtPhaseIbs phaseIbs) {
        this.phaseData = phaseIbs.phaseData();
        this.burnin = phaseData.it() < phaseData.fpd().par().burnin();
        this.estPhase = phaseData.estPhase();
        this.markers = phaseData.fpd().stage1TargGT().markers();
        this.nMarkers = markers.size();
        this.maxStates = phaseData.fpd().par().phase_states();
        this.states = new BasicPhaseStates(phaseIbs, maxStates);
        this.lrList = new FloatList(200);

        this.refAlleles = new ArrayList<>();
        this.mismatches = new byte[3][nMarkers][maxStates];
        this.pMismatch = phaseData.pMismatch();
        this.emProbs = new float[] {1.0f - pMismatch, pMismatch};

        this.fwd = new float[3][maxStates];
        this.bwd = new float[3][maxStates];
        this.sum = new float[3];
        this.missProbs1 = new ArrayList<>();
        this.missProbs2 = new ArrayList<>();
        this.bwdHet1 = new ArrayList<>();
        this.bwdHet2 = new ArrayList<>();
    }

    /**
     * Returns the number of target samples.
     * @return the number of target samples
     */
    public int nTargSamples() {
        return phaseData.fpd().targGT().nSamples();
    }

    /**
     * Estimates and stores the phased haplotypes for the specified sample
     * @param sample a sample index
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nTargSamples()}
     */
    public void phase(int sample) {
        SamplePhase phase = estPhase.get(sample);
        swapHaps = false;
        swapCnt = 0;
        int nUnph = phase.unphased().size();
        int nMiss = phase.missing().size();
        if (nMiss>0 || nUnph>0) {
            lrList.clear();
            ensureCapacity(nUnph, nMiss);
            MarkerCluster mc = new MarkerCluster(phaseData, sample);
            missIndex = mc.nMissingGTClusters();
            nStates = states.ibsStates(sample, mc, refAlleles, mismatches);
            bwdAlg(mc);
            fwdAlg(phase, mc);
            updatePhase(sample, phase);
            estPhase.set(sample, phase);
        }
        nUnphHets.addAndGet(nUnph);
        nSwaps.addAndGet(swapCnt);
    }

    private void ensureCapacity(int nUnph, int nMiss) {
        if (refAlleles.size()<nMiss) {
            for (int j=refAlleles.size(); j<nMiss; ++j) {
                refAlleles.add(new int[maxStates]);
                missProbs1.add(new float[maxStates]);
                missProbs2.add(new float[maxStates]);
            }
        }
        if (bwdHet1.size()<nUnph) {
            for (int j=bwdHet1.size(); j<nUnph; ++j) {
                bwdHet1.add(new float[maxStates]);
                bwdHet2.add(new float[maxStates]);
            }
        }
    }

    private void fwdAlg(SamplePhase phase, MarkerCluster mc) {
        IntArray unph = mc.unphClusters();
        Arrays.fill(fwd[0], 0, nStates, 1.0f/nStates);
        sum[0] = 1.0f;
        int startClust = 0;
        for (int j=0, n=unph.size(); j<n; ++j) {
            int endClust = unph.get(j);
            fwdAlg(phase, mc, startClust, endClust);
            phaseHet(j);
            startClust = endClust;
        }
        if (startClust<mc.nClusters()) {
            fwdAlg(phase, mc, startClust, mc.nClusters());
        }
    }

    private void fwdAlg(SamplePhase phase, MarkerCluster mc, int startClust,
            int endClust) {
        if (swapHaps) {
            swapHaps(phase, mc, startClust, endClust);
        }
        System.arraycopy(fwd[0], 0, fwd[1], 0, nStates);
        System.arraycopy(fwd[0], 0, fwd[2], 0, nStates);
        sum[1] = sum[2] = sum[0];
        FloatArray pClustRecomb = mc.pRecomb();
        for (int c=startClust; c<endClust; ++c)  {
            float pRec = pClustRecomb.get(c);
            emProbs[1] = (mc.clusterEnd(c) - mc.clusterStart(c))*pMismatch;
            emProbs[0] = 1.0f - emProbs[1];
            sum[0] = HmmUpdater.fwdUpdate(fwd[0], sum[0], pRec, emProbs, mismatches[0][c], nStates);
            sum[1] = HmmUpdater.fwdUpdate(fwd[1], sum[1], pRec, emProbs, mismatches[1][c], nStates);
            sum[2] = HmmUpdater.fwdUpdate(fwd[2], sum[2], pRec, emProbs, mismatches[2][c], nStates);
            if (mc.clustHasMissingGT(c)) {
                imputeAlleles(phase, mc, c);
            }
        }
    }

    private void swapHaps(SamplePhase phase, MarkerCluster mc, int startClust, int endClust) {
        for (int c=startClust; c<endClust; ++c) {
            byte[] tmpMatch = mismatches[1][c];
            mismatches[1][c] = mismatches[2][c];
            mismatches[2][c] = tmpMatch;
        }
        phase.swapHaps(mc.clusterStart(startClust), mc.clusterEnd(endClust-1));
    }

    private void imputeAlleles(SamplePhase phase, MarkerCluster mc, int cluster) {
        if (mc.clustHasMissingGT(cluster)) {
            if (swapHaps) {
                float[] tmp = missProbs1.get(missIndex);
                missProbs1.set(missIndex, missProbs2.get(missIndex));
                missProbs2.set(missIndex, tmp);
            }
            float[] stateProbs1 = missProbs1.get(missIndex);
            float[] stateProbs2 = missProbs2.get(missIndex);
            int[] refAl = refAlleles.get(missIndex);
            for (int k=0; k<nStates; ++k) {
                stateProbs1[k] *= fwd[1][k];
                stateProbs2[k] *= fwd[2][k];
            }

            assert (mc.clusterEnd(cluster) - mc.clusterStart(cluster))==1;
            int marker = mc.clusterStart(cluster);
            imputeAlleles(phase, marker, stateProbs1, stateProbs2, refAl);
            ++missIndex;
        }
    }

    private void imputeAlleles(SamplePhase phase, int marker,
            float[] stateProbs1, float[] stateProbs2, int[] refAl) {
        int nAlleles = markers.marker(marker).nAlleles();
        float[] alFreq1 = new float[nAlleles];
        float[] alFreq2 = new float[nAlleles];
        for (int k=0; k<nStates; ++k) {
            alFreq1[refAl[k]] += stateProbs1[k];
            alFreq2[refAl[k]] += stateProbs2[k];
        }
        int a1 = 0;
        int a2 = 0;
        for (int j=1; j<nAlleles; ++j) {
            if (alFreq1[j]>alFreq1[a1]) {
                a1 = j;
            }
            if (alFreq2[j]>alFreq2[a2]) {
                a2 = j;
            }
        }
        phase.setAllele1(marker, a1);
        phase.setAllele2(marker, a2);
    }

    private void bwdAlg(MarkerCluster mc) {
        IntArray unph = mc.unphClusters();
        Arrays.fill(bwd[0], 0, nStates, 1.0f/nStates);
        int endClust = mc.nClusters() - 1;
        if (mc.clustHasMissingGT(endClust)) {
            --missIndex;
            System.arraycopy(bwd[0], 0, missProbs1.get(missIndex), 0, nStates);
            System.arraycopy(bwd[0], 0, missProbs2.get(missIndex), 0, nStates);
        }
        for (int j=unph.size()-1; j>=0; --j) {
            int startClust = unph.get(j)-1;
            assert startClust>=0; // first het is not in unphased list
            bwdAlg(mc, startClust, endClust);
            System.arraycopy(bwd[1], 0, bwdHet1.get(j), 0, nStates);
            System.arraycopy(bwd[2], 0, bwdHet2.get(j), 0, nStates);
            endClust = startClust;
        }
        bwdAlg(mc, 0, endClust);
    }

    private void bwdAlg(MarkerCluster mc, int startClust, int endClust) {
        FloatArray pClustRecomb = mc.pRecomb();
        System.arraycopy(bwd[0], 0, bwd[1], 0, nStates);
        System.arraycopy(bwd[0], 0, bwd[2], 0, nStates);
        for (int c=endClust-1; c>=startClust; --c) {
            int cP1 = c + 1;
            float pRec = pClustRecomb.get(cP1);
            emProbs[1] = (mc.clusterEnd(c) - mc.clusterStart(c))*pMismatch;
            emProbs[0] = 1f - emProbs[1];
            HmmUpdater.bwdUpdate(bwd[0], pRec, emProbs, mismatches[0][cP1], nStates);
            HmmUpdater.bwdUpdate(bwd[1], pRec, emProbs, mismatches[1][cP1], nStates);
            HmmUpdater.bwdUpdate(bwd[2], pRec, emProbs, mismatches[2][cP1], nStates);
            if (mc.clustHasMissingGT(c)) {
                --missIndex;
                System.arraycopy(bwd[1], 0, missProbs1.get(missIndex), 0, nStates);
                System.arraycopy(bwd[2], 0, missProbs2.get(missIndex), 0, nStates);
            }
        }
    }

    private void phaseHet(int hetIndex) {
        float[] b1 = bwdHet1.get(hetIndex);
        float[] b2 = bwdHet2.get(hetIndex);
        float p11 = 0.0f;
        float p12 = 0.0f;
        float p21 = 0.0f;
        float p22 = 0.0f;
        for (int k=0; k<nStates; ++k) {
            p11 += fwd[1][k]*b1[k];
            p12 += fwd[1][k]*b2[k];
            p21 += fwd[2][k]*b1[k];
            p22 += fwd[2][k]*b2[k];
        }
        boolean lastSwapHaps = swapHaps;
        float num = (p11*p22);
        float den = (p12*p21);
        swapHaps = num < den;
        if (swapHaps!=lastSwapHaps) {
            ++swapCnt;
        }
        lrList.add(swapHaps ? den/num : num/den);
    }

    private void updatePhase(int sample, SamplePhase phase) {
        IntArray prevUnph = phase.unphased();
        if (prevUnph.size()>0 && burnin==false) {
            float leaveUnphasedProp = phaseData.leaveUnphasedProp(sample);
            IntList nextUnph = new IntList();
            float threshold = threshold(lrList, leaveUnphasedProp);
            for (int j=0, n=prevUnph.size(); j<n; ++j) {
                if (lrList.get(j) < threshold) {
                    nextUnph.add(prevUnph.get(j));
                }
            }
            phase.setUnphased(IntArray.packedCreate(nextUnph, nMarkers));
        }
    }

    private static float threshold(FloatList lrList, float leaveUnphasedProp) {
        float[] lra = lrList.toArray();
        Arrays.sort(lra);
        int rank = (int) Math.floor(leaveUnphasedProp*lra.length + 0.5f);
        return lra[rank<lra.length ? rank : (lra.length-1)];
    }
}
