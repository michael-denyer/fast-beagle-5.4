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
import ints.IntArray;
import ints.IntList;
import ints.WrappedIntArray;

/**
 * <p>Class {@code MarkerCluster} represents a partition of markers into
 * contiguous marker clusters.</p>
 *
 * <p>Instances of class {@code MarkerCluster} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class MarkerCluster {

    private final int[] clusterToEnd;
    private final IntArray hetClusters;
    private final int nMissingGTClusters;
    private final boolean[] clustHasMissingGT;
    private final FloatArray pRecomb;

    /**
     * Constructs a new {@code MarkerCluster} instance from the specified data.
     *
     * @param phaseData the input data for the next genotype phasing iteration
     * @param sample a sample index
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= phaseData.targGT().nSamples()}
     * @throws NullPointerException if {@code phaseData == null}
     */
    public MarkerCluster(PhaseData phaseData, int sample) {
        SamplePhase samplePhase = phaseData.estPhase().get(sample);
        this.clusterToEnd = samplePhase.clustEnds();
        this.hetClusters = unphHetClusters(samplePhase, clusterToEnd);
        boolean[] hasMissingGT = new boolean[clusterToEnd.length];
        this.nMissingGTClusters = setClustHasMissingGT(samplePhase.missing(),
                clusterToEnd, hasMissingGT);
        this.clustHasMissingGT = hasMissingGT;
        this.pRecomb = pClustRecomb(phaseData.pRecomb(), clusterToEnd);
    }

    private static IntArray unphHetClusters(SamplePhase samplePhase,
            int[] clustToEnd) {
        IntArray unph = samplePhase.unphased();
        int nUnph = unph.size();
        IntList hetClusters = new IntList(nUnph);
        int unphIndex = 0;
        int nextUnph = unphIndex<nUnph ? unph.get(unphIndex++) : Integer.MAX_VALUE;
        for (int j=0; j<clustToEnd.length; ++j) {
            int clustEnd = clustToEnd[j];
            if (nextUnph<clustEnd) {
                hetClusters.add(j);
                nextUnph = unphIndex<nUnph ? unph.get(unphIndex++) : Integer.MAX_VALUE;
                while (nextUnph<clustEnd) {
                    nextUnph = unphIndex<nUnph ? unph.get(unphIndex++) : Integer.MAX_VALUE;
                }
            }
        }
        return new WrappedIntArray(hetClusters);
    }

    private static int setClustHasMissingGT(IntArray missGTList,
            int[] cluster2End, boolean[] clustHasMissingGT) {
        int missCnt = 0;
        int i = 0;
        int nMissGT = missGTList.size();
        for (int c=0; c<cluster2End.length; ++c) {
            int end = cluster2End[c];
            for ( ; i<nMissGT && missGTList.get(i)<end; ++i) {
                clustHasMissingGT[c] = true;
            }
            if (clustHasMissingGT[c]) {
                ++missCnt;
            }
        }
        return missCnt;
    }

    private static FloatArray pClustRecomb(FloatArray pRecomb,
            int[] cluster2End) {
        int nClusters = cluster2End.length;
        float[] pClustRecomb = new float[nClusters];
        int start = 0;
        for (int j=1; j<nClusters; ++j) {
            int end = cluster2End[j];
            float pNoRecomb = 1.0f;
            for (int k=start; k<end; ++k) {
                pNoRecomb *= (1.0f - pRecomb.get(k));
            }
            pClustRecomb[j] = 1.0f - pNoRecomb;
            start = end;
        }
        return new FloatArray(pClustRecomb);
    }

    /**
     * Returns the number of clusters
     * @return the number of clusters
     */
    public int nClusters() {
        return clusterToEnd.length;
    }

    /**
     * Returns the inclusive start marker for the cluster.
     * @param index a cluster index
     * @return the inclusive start marker for the cluster
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nClusteres()}
     */
    public int clusterStart(int index) {
        return index==0 ? 0 : clusterToEnd[index-1];
    }

    /**
     * Returns the exclusive end marker for the cluster.
     * @param index a cluster index
     * @return the exclusive marker for the cluster
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nClusteres()}
     */
    public int clusterEnd(int index) {
        return clusterToEnd[index];
    }

    /**
     * Return a {@code FloatArray} of size {@code this.nClusters()}
     * whose {@code k}-th element is the probability of transitioning to a
     * random HMM state between the {@code k}-th cluster and the
     * previous cluster.
     * @return a {@code FloatArray} of size {@code this.nClusters()}
     * whose {@code k}-th element is the probability of transitioning to a
     * random HMM state between the {@code k}-th cluster and the
     * previous cluster
     */
    public FloatArray pRecomb() {
        return pRecomb;
    }

    /**
     * Returns a sorted list of cluster indices in increasing order for which
     * the cluster contains an unphased heterozygote.
     * @return a sorted list of cluster indices in increasing order for which
     * the cluster contains an unphased heterozygote
     */
    public IntArray unphClusters() {
        return hetClusters;
    }

    /**
     * Returns {@code true} if the cluster has at least one missing genotype,
     * and returns {@code false} otherwise.
     * @param index a cluster index
     * @return {@code true} if the cluster has at least one missing genotype
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nClusters()}
     */
    public boolean clustHasMissingGT(int index) {
        return clustHasMissingGT[index];
    }

    /**
     * Returns the number of clusters containing at least one missing genotype.
     * @return the number of clusters with at least one missing genotype
     */
    public int nMissingGTClusters() {
        return nMissingGTClusters;
    }
}
