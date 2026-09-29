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

import vcf.Steps;
import blbutil.Utilities;
import ints.IntIntMap;
import ints.IntList;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.stream.IntStream;
import beagleutil.CompHapSegment;
import vcf.XRefGT;

/**
 * <p>Class {@code LowFreqPhaseStates} has methods for constructing a Li and
 * Stephens HMM for a target haplotype.  The resulting HMM states are
 * enriched for reference haplotypes carrying low frequency variants.
 * </p>
 * <p>Instances of {@code LowFreqPhaseStates} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class LowFreqPhaseStates {

    private static final int NIL = -103;
    private final LowFreqPhaseIbs ibsHaps;
    private final PhaseData phaseData;
    private final Steps steps;
    private final XRefGT allHaps;
    private final int nMarkers;
    private final int maxStates;
    private final int minSteps;

    private final IntIntMap hapToEnd;
    private final PriorityQueue<CompHapSegment> q;

    private final IntList[] compositeHapToHap;
    private final IntList[] compositeHapToEnd;
    private final int[] segmentIndex;
    private final int[] compHapToHap;
    private final int[] compHapToEnd;

    /**
     * Constructs a new {@code LowFreqPhaseStates} object from the specified
     * data.
     * @param ibsHaps the IBS haplotype segments
     * @param maxStates the maximum number of composite reference
     * haplotypes that will be constructed
     * @throws IllegalArgumentException if {@code maxStates < 1}
     * @throws NullPointerException if {@code ibsHaps == null}
     */
    public LowFreqPhaseStates(LowFreqPhaseIbs ibsHaps, int maxStates) {
        if (maxStates < 1) {
            throw new IllegalArgumentException(String.valueOf(maxStates));
        }
        this.ibsHaps = ibsHaps;
        this.phaseData = ibsHaps.phaseData();
        this.steps = phaseData.fpd().stage1Steps();
        this.allHaps = ibsHaps.allHaps();
        this.nMarkers = allHaps.nMarkers();
        this.maxStates = maxStates;
        float phaseStep = phaseData.fpd().ibsStep();
        this.minSteps = Math.max(200, (int) Math.ceil(1.0f/phaseStep)); // 200 steps and 1 cM
        this.hapToEnd = new IntIntMap(maxStates);
        this.q = new PriorityQueue<>(maxStates);

        this.compositeHapToHap = IntStream.range(0, maxStates)
                .mapToObj(j -> new IntList())
                .toArray(IntList[]::new);
        this.compositeHapToEnd = IntStream.range(0, maxStates)
                .mapToObj(j -> new IntList())
                .toArray(IntList[]::new);
        this.segmentIndex = new int[maxStates];
        this.compHapToHap = new int[maxStates];
        this.compHapToEnd = new int[maxStates];
    }

    /**
     * Returns the number of target haplotypes.
     * @return the number of target haplotypes
     */
    public int nTargHaps() {
        return phaseData.fpd().targGT().nHaps();
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return phaseData.fpd().targGT().nMarkers();
    }

    /**
     * Returns the maximum number of HMM states at a marker.
     * @return the maximum number of HMM states at a marker
     */
    public int maxStates() {
        return maxStates;
    }

    /**
     * Stores the Li and Stephens HMM for the specified target
     * haplotype in the specified arrays.  The haplotype for the
     * {@code j}-th state at the {@code m}-th marker is stored
     * in {@code haps[m][j]}.  The number of allele mismatches (0 or 1)
     * between the haplotype for the {@code j}-th state and the
     * target haplotype at the {@code m}-th marker is stored in
     * {@code nMismatches[m][j]}.
     * The number of HMM states states at each marker is returned.
     * @param targHap the haplotype index
     * @param haps the two-dimensional array in which the
     * haplotype for each HMM state will be stored
     * @param nMismatches the two-dimensional array in which the number
     * of allele mismatches (0 or 1) for each HMM state will be stored
     * @return the number of HMM states at each marker
     *
     * @throws IndexOutOfBoundsException if
     * {@code targHap < 0 || targHap >= this.nTargHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code haps.length < this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code nMismatches.length < this.nMarkers()}
     * @throws IndexOutOfBoundsException if {@code haps[m].length}
     * is less than the number of HMM states for any marker {@code m}
     * satisfying {@code (0 <= m && m < haps.length)}
     * @throws IndexOutOfBoundsException if {@code nMismatches[m].length}
     * is less than the number of HMM states for any marker {@code m}
     * satisfying {@code (0 <= m && m < nMismatches.length)}
     * @throws NullPointerException if any array is {@code null}
     */
    public int ibsStates(int targHap, int[][] haps, byte[][] nMismatches) {
        int nCompHaps = setCompRefHaps(targHap);
        copyData(targHap, nCompHaps, haps, nMismatches);
        return nCompHaps;
    }

    private int setCompRefHaps(int targHap) {
        q.clear();
        hapToEnd.clear();
        for (int j=0, n=maxStates; j<n; ++j) {
            compositeHapToHap[j].clear();
            compositeHapToEnd[j].clear();
        }
        for (int step=0, n=steps.size(); step<n; ++step) {
            //ibsHaps.ibsHaps(targHap, step, ibsHapList);
            addIbsHap(ibsHaps.fwdIbsHap(targHap, step), step);
            addIbsHap(ibsHaps.bwdIbsHap(targHap, step), step);
        }
        if (q.isEmpty()) {
            fillQWithRandomHaps(targHap);
        }
        int nCompHaps = setFinalRefSegs();
        return nCompHaps;
   }

    private void addIbsHap(int ibsHap, int step) {
        if (ibsHap<0) {
            return;
        }
        if (hapToEnd.get(ibsHap, NIL)==NIL) { // hap is not currently in q
            updateHeadOfQ();
            if (q.size()==maxStates
                    || (q.isEmpty()==false && (step - q.peek().ibsStep()) >= minSteps)) {
                CompHapSegment head = q.poll();
                int index = head.compHapIndex();
                int prevHap = head.hap();
                int nextStart = steps.start((head.ibsStep() + step) >>> 1);
                hapToEnd.remove(prevHap);
                compositeHapToHap[index].add(ibsHap);      // hap of new segment
                compositeHapToEnd[index].add(nextStart);   // end of old segment

                head.updateSegment(ibsHap, nextStart, step);
                q.add(head);
            }
            else {
                int index = q.size();
                compositeHapToHap[index].add(ibsHap);            // hap of new segment
                q.add(new CompHapSegment(ibsHap, 0, step, index));
            }
        }
        hapToEnd.put(ibsHap, step);
    }

    private void updateHeadOfQ() {
        CompHapSegment head = q.peek();
        if (head!=null) {
            int latestEnd = hapToEnd.get(head.hap(), NIL);
            while (head.ibsStep()!=latestEnd) {
                head = q.poll();
                head.updateStep(latestEnd);
                q.offer(head);
                head = q.peek();
                latestEnd = hapToEnd.get(head.hap(), NIL);
            }
        }
    }

    private int setFinalRefSegs() {
        int nCompHaps = q.size();
        CompHapSegment head = q.poll();
        while (head!=null) {
            int compHap = head.compHapIndex();
            compositeHapToEnd[compHap].add(nMarkers); // add missing end of last segment
            segmentIndex[compHap] = 0;
            compHapToHap[compHap] = compositeHapToHap[compHap].get(0);
            compHapToEnd[compHap] = compositeHapToEnd[compHap].get(0);
            head = q.poll();
        }
        return nCompHaps;
    }

    private void copyData(int targHap, int nCompHaps, int[][] haps, byte[][] nMismatches) {
        for (int m=0; m<nMarkers; ++m) {
            int obsAllele = allHaps.allele(m, targHap);
            for (int j=0; j<nCompHaps; ++j) {
                if (m==compHapToEnd[j]) {
                    ++segmentIndex[j];
                    compHapToHap[j] = compositeHapToHap[j].get(segmentIndex[j]);
                    compHapToEnd[j] = compositeHapToEnd[j].get(segmentIndex[j]);
                }
                int refHap = compHapToHap[j];
                haps[m][j] = refHap;
                nMismatches[m][j] = allHaps.allele(m, refHap)==obsAllele
                        ? (byte) 0 : (byte) 1;
            }
        }
    }

    private void fillQWithRandomHaps(int hap) {
        assert q.isEmpty();
        int nHaps = allHaps.nHaps();
        int nStates = Math.min(nHaps-2, maxStates);
        if (nStates<=0) {
            Utilities.exit("ERROR: there is only one sample");
        }
        else {
            Random rand = new Random(phaseData.seed() + hap);
            int step = steps.size()-1;
            int sample = hap>>1;
            for (int i=0; i<nStates; ++i) {
                int h = rand.nextInt(nHaps);
                while ((h>>1)==sample) {
                    h = rand.nextInt(nHaps);
                }
                compositeHapToHap[q.size()].add(h);
                q.add(new CompHapSegment(h, 0, step, i));
            }
        }
    }
}
