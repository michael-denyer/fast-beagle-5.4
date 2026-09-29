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

import blbutil.BitArray;
import blbutil.DoubleArray;
import ints.IntArray;
import ints.IntList;
import java.util.Arrays;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import vcf.Markers;

/**
 * <p>Each instance of class {@code SamplePhase} stores an estimated haplotype
 * pair for a sample, the list of markers with missing genotypes for the sample,
 * a list of markers whose genotype phase with respect to the preceding
 * heterozygote genotype is considered to be uncertain for the sample, and
 * a set of marker clusters for the sample.
 * </p>
 * <p>Instances of class {@code SamplePhase} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class SamplePhase {

    private final Markers markers;
    private BitArray hap1;
    private BitArray hap2;
    private IntArray unphased;
    private final IntArray missing;
    private final byte[] clustSize;

    /**
     * Constructs a new {@code SamplePhase} instance from the specified data.
     * @param markers the list of markers
     * @param genPos the genetic positions of the specifed markers
     * @param hap1 the list of alleles on the first haplotype
     * @param hap2 the list of alleles on the second haplotype
     * @param unphased the indices of markers whose genotype phase with respect
     * to the preceding heterozygote is unknown
     * @param missing the indices of markers whose genotype is missing
     * @throws IllegalArgumentException if
     * {@code genPos.size() != markers.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code hap1.length != markers.nMarkers()
     * || hap2.length != markers.nMarkers()}
     * @throws IllegalArgumentException if the specified {@code unphased} or
     * {@code missing} list is not a strictly increasing list of
     * marker indices between 0 (inclusive) and {@code markers.nMarkers()}
     * (exclusive)
     * @throws NullPointerException if any argument is {@code null}
     */
    public SamplePhase(Markers markers, DoubleArray genPos,
            int[] hap1, int[] hap2, IntArray unphased, IntArray missing) {
        int nMarkers = markers.size();
        if (nMarkers!=genPos.size()) {
            throw new IllegalArgumentException(String.valueOf(genPos.size()));
        }
        if (hap1.length!=nMarkers) {
            throw new IllegalArgumentException(String.valueOf(hap1.length));
        }
        if (hap2.length!=nMarkers) {
            throw new IllegalArgumentException(String.valueOf(hap2.length));
        }
        checkIncreasing(unphased, nMarkers);
        checkIncreasing(missing, nMarkers);
        this.markers = markers;
        this.hap1 = new BitArray(markers.sumHapBits());
        this.hap2 = new BitArray(markers.sumHapBits());
        markers.allelesToBits(hap1, this.hap1);
        markers.allelesToBits(hap2, this.hap2);
        this.unphased = unphased;
        this.missing = missing;
        float maxClusterCM = 0.005f;
        this.clustSize = clustSize(hap1, hap2, missing, genPos, maxClusterCM);
    }

    private static void checkIncreasing(IntArray ia, int nMarkers) {
        int last = -1;
        for (int j=0, n=ia.size(); j<n; ++j) {
            if (ia.get(j)<=last) {
                throw new IllegalArgumentException(ia.toString());
            }
            last = ia.get(j);
        }
        if (last>=nMarkers) {
            throw new IllegalArgumentException(ia.toString());
        }
    }

    private static byte[] clustSize(int[] hap1, int[] hap2, IntArray missing,
            DoubleArray genPos, float maxCM)  {
        IntList clustSizes = new IntList(1<<12);
        int nMarkers = genPos.size();
        double maxClustEnd = genPos.get(0) + maxCM;
        boolean prevIsMissOrHet = false;
        int lastEnd = 0;
        int missIndex = 0;
        int nextMiss = missIndex<missing.size() ? missing.get(missIndex++) : -1;
        for (int m=0; m<nMarkers; ++m) {
            int size = m - lastEnd;
            boolean isMissing = m==nextMiss;
            if (isMissing) {
                nextMiss = missIndex<missing.size() ? missing.get(missIndex++) : -1;
            }
            boolean isMissOrHet = isMissing || hap1[m]!=hap2[m];
            if (prevIsMissOrHet || isMissOrHet || genPos.get(m)>maxClustEnd || size==255) {
                if (m>0) {
                    clustSizes.add(size);
                    maxClustEnd = genPos.get(m) + maxCM;
                    lastEnd = m;
                }
            }
            prevIsMissOrHet = isMissOrHet;
        }
        clustSizes.add(nMarkers - lastEnd);
        return toByteArray(clustSizes);
    }

    private static byte[] toByteArray(IntList intList) {
        byte[] ba = new byte[intList.size()];
        for (int j=0; j<ba.length; ++j) {
            ba[j] = (byte) intList.get(j);
        }
        return ba;
    }

    /**
     * Returns the (exclusive) end marker indices of each marker cluster.
     * The returned list is sorted in increasing order.
     * @return the (exclusive) end marker indices of each marker cluster
     */
    public int[] clustEnds() {
        int[] clustEnds = new int[clustSize.length];
        int cumSum = 0;
        for (int j=0; j<clustSize.length; ++j) {
            cumSum += (clustSize[j] & 0xff); // convert unsigned byte to integer
            clustEnds[j] = cumSum;
        }
        return clustEnds;
    }

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    public Markers markers() {
        return markers;
    }

    /**
     * Returns a list of marker indices in increasing order for which
     * the genotype is missing.
     * @return a list of marker indices in increasing order for which
     * the genotype is missing
     */
    public IntArray missing() {
        return missing;
    }

    /**
     * Returns a list of marker indices in increasing order whose genotype
     * phase with respect to the preceding non-missing heterozygote genotype
     * is unknown.
     * @return a list of markers indices in increasing order whose genotype
     * phase with respect to the preceding non-missing heterozygote genotype
     * is unknown
     */
    public IntArray unphased() {
        return unphased;
    }

    /**
     * Sets the list of markers whose genotype phase with respect to
     * the preceding non-missing heterozygote genotype is unknown.
     * @param unphased a list of markers whose genotype phase with respect to
     * the preceding non-missing heterozygote genotype is unknown
     * @throws IllegalArgumentException if the specified list or marker
     * indices is not a strictly increasing list of indices between 0
     * (inclusive) and {@code this.markers().nMarkers()} (exclusive)
     * @throws NullPointerException if {@code unphased == null}
     */
    public void setUnphased(IntArray unphased) {
        checkIncreasing(unphased, markers.size());
        this.unphased = unphased;
    }

    /**
     * Copies the stored haplotypes to the specified {@code BitList} objects
     * @param hap1 a {@code BitList} in which the sample's first haplotype's
     * alleles will be  stored
     * @param hap2 a {@code BitList} in which the sample's second haplotype's
     * alleles will be  stored
     * @throws IllegalArgumentException if
     * {@code hap1.size() != this.markers().sumHaplotypeBits()}
     * @throws IllegalArgumentException if
     * {@code hap2.size()!= this.markers().sumHaplotypeBits()}
     * @throws NullPointerException if {@code hap1 == null || hap2 == null}
     */
    public void getHaps(BitArray hap1, BitArray hap2) {
        int nBits = markers.sumHapBits();
        if (hap1.size() != nBits || hap2.size() != nBits) {
            throw new IllegalArgumentException("inconsistent data");
        }
        hap1.copyFrom(this.hap1, 0, this.hap1.size());
        hap2.copyFrom(this.hap2, 0, this.hap2.size());
    }


    /**
     * Returns the allele on the first haplotype for the specified marker.
     * @param marker the marker index
     * @return the allele on the first haplotype for the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.markers().nMarkers()}
     */
    public int allele1(int marker) {
       return markers.allele(hap1, marker);
    }

    /**
     * Returns the allele on the second haplotype for the specified marker.
     * @param marker the marker index
     * @return the allele on the second haplotype for the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.markers().nMarkers()}
     */
    public int allele2(int marker) {
        return markers.allele(hap2, marker);
    }

    /**
     * Sets the allele on the first haplotype for the specified marker
     * to the specified allele
     * @param marker the marker index
     * @param allele the allele
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.markers().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.markers().marker(marker).nAlleles()}
     */
    public void setAllele1(int marker, int allele) {
        markers.setAllele(marker, allele, hap1);
    }

    /**
     * Sets the allele on the second haplotype for the specified marker
     * to the specified allele
     * @param marker the marker index
     * @param allele the allele
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.markers().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.markers().marker(marker).nAlleles()}
     */
    public void setAllele2(int marker, int allele) {
        markers.setAllele(marker, allele, hap2);
    }

    /**
     * Swaps the alleles of the two haplotypes in the specified range of
     * markers.
     * @param start the start marker index (inclusive)
     * @param end the end marker index (exclusive)
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start > end || start >= this.markers().nMarkers()}
     */
    public void swapHaps(int start, int end) {
        int startBit = markers.sumHapBits(start);
        int endBit = markers.sumHapBits(end);
        BitArray.swapBits(hap1, hap2, startBit, endBit);
    }

    /**
     * Returns the first haplotype.  The haplotype is encoded with the
     * {@code this.markers().allelesToBits()} method.
     * @return the first haplotype
     */
    public BitArray hap1() {
        return new BitArray(this.hap1);
    }

    /**
     * Returns the second haplotype.  The haplotype is encoded with the
     * {@code this.markers().allelesToBits()} method.
     * @return the second haplotype
     */
    public BitArray hap2() {
        return new BitArray(this.hap2);
    }

    /**
     * Returns the current estimated phased genotypes.  This method converts
     * column-major data into row-major data.
     * @param estPhase the current estimated phased genotypes for each target
     * sample
     * @return the current estimated phased genotypes
     * @throws NullPointerException if {@code estPhase == null}
     */
    public static BitArray[] toBitLists(EstPhase estPhase) {
        int nThreads = estPhase.fpd().par().nthreads();
        int nMarkers = estPhase.fpd().stage1TargGT().nMarkers();
        int nRecsPerBatch = (nMarkers + nThreads - 1)/nThreads;
        while (nRecsPerBatch>4096) {
            nRecsPerBatch = (nRecsPerBatch+1) >> 1;
        }
        int stepSize = nRecsPerBatch;
        int nSteps = (nMarkers + (stepSize-1)) / stepSize;
        return IntStream.range(0, nSteps)
                .parallel()
                .boxed()
                .flatMap(step -> bitLists(estPhase, step, stepSize))
                .toArray(BitArray[]::new);
    }

    private static Stream<BitArray> bitLists(EstPhase estPhase, int step, int stepSize) {
        int nSamples = estPhase.fpd().targGT().nSamples();
        int nHaps = nSamples<<1;
        Markers markers = estPhase.fpd().stage1TargGT().markers();
        int mStart = step*stepSize;
        int mEnd = Math.min(mStart + stepSize, markers.size());
        BitArray[] bitLists = IntStream.range(mStart, mEnd)
                .mapToObj(j -> new BitArray(nHaps*markers.marker(j).bitsPerAllele()))
                .toArray(BitArray[]::new);
        int[] bitsPerAllele = IntStream.range(mStart, mEnd)
                .map(m -> markers.marker(m).bitsPerAllele())
                .toArray();
        for (int s=0; s<nSamples; ++s) {
            SamplePhase sampPhase = estPhase.get(s);
            int h1 = s<<1;
            int h2 = h1 | 0b1;
            int inBit1 = markers.sumHapBits(mStart);
            int inBit2 = markers.sumHapBits(mStart);
            for (int m=mStart; m<mEnd; ++m) {
                int mOffset = m - mStart;
                int nBits = bitsPerAllele[mOffset];
                int startOutBit1 = h1*nBits;
                int startOutBit2 = h2*nBits;
                for (int i=0; i<nBits; ++i) {
                    if (sampPhase.hap1.get(inBit1++)) {
                        bitLists[mOffset].set(startOutBit1 + i);
                    }
                    if (sampPhase.hap2.get(inBit2++)) {
                        bitLists[mOffset].set(startOutBit2 + i);
                    }
                }
            }
        }
        return Arrays.stream(bitLists);
    }
}
