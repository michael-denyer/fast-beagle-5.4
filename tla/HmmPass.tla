------------------------------- MODULE HmmPass -------------------------------
(***************************************************************************)
(* One phase_baum1 (pb) phasing S samples in turn, from                    *)
(* src/phase/phase_baum1.c: begin_sample, bwd_begin, bwd_step with         *)
(* bwd_end_step, fwd_begin, fwd_step with phase_het and swap_haps,         *)
(* fwd_end_step with impute_alleles, and end_sample with update_phase and  *)
(* threshold. Each sample has clusters 0 .. N-1, each missing-genotype     *)
(* cluster ("MISS"), homozygous, phased het ("HET") or unphased het        *)
(* ("UNPH"): marker_cluster.c gives each unphased marker a cluster of its  *)
(* own (ClusterLifecycle checks it), so the UNPH clusters stand for        *)
(* sp->unphased. Init chooses pb->burnin, which phase_baum1_init sets for  *)
(* the whole run, and the cluster types, with a phased het before every    *)
(* UNPH, since PbwtPhaser.indices (pbwt_phaser.c) never lists a sample's   *)
(* first heterozygote. At each unphased het the environment chooses num <  *)
(* den and the likelihood ratio phase_het stores: 1, 2, or NaN (0 / 0,     *)
(* only when num < den is false); at end_sample, the rank                  *)
(* floor(leave_unph_prop * n + 0.5), from 0 to n because the proportion is *)
(* at most 1. A read of lr_sorted outside [0, n) returns any value. Ratios *)
(* sort as Float.compare does, NaN last. The backward values are not       *)
(* modeled, only the indices the backward pass stores them under and the   *)
(* forward pass reads them from. Each step is one bwd_step with            *)
(* bwd_end_step, one fwd_step with fwd_end_step, or a sample's begin or    *)
(* end.                                                                    *)
(***************************************************************************)
EXTENDS Integers, Sequences, FiniteSets, TLC

CONSTANTS N, S

ASSUME N \in Nat \ {0, 1} /\ S \in Nat \ {0}

VARIABLES
    burnin,     \* pb->burnin
    types0,     \* each sample's cluster types when the run started
    types,      \* each sample's cluster types; UNPH becomes HET in update_phase
    snap,       \* the current sample's types at begin_sample
    pc,         \* "begin", "bwd", "fwd", "end" or "done"
    s,          \* the current sample
    ran,        \* the samples whose begin_sample returned true
    swapHaps,   \* pb->swap_haps
    swapCnt,    \* pb->swap_cnt
    lr,         \* pb->lr[0 .. n_lr)
    outs,       \* per sample, each phase_het's swap_haps
    missIdx,    \* pb->miss_index, shared by both passes
    bc, bUnph,  \* bwd_pass c and unph_index
    fc, fUnph,  \* fwd_pass c and unph_index
    hetW, missW,  \* per cluster, the index the backward pass stored it under
    hetR, missR,  \* per cluster, the index the forward pass read it from
    swaps,      \* per cluster, swap_haps calls covering it in this sample
    rateSwaps, rateUnph,  \* swap_rate n_swaps and n_unph_hets
    upd         \* the last update_phase: sample, types before, index read, threshold

vars == <<burnin, types0, types, snap, pc, s, ran, swapHaps, swapCnt, lr, outs, missIdx,
          bc, bUnph, fc, fUnph, hetW, missW, hetR, missR, swaps, rateSwaps, rateUnph, upd>>

Cl == 0 .. N - 1
Samples == 1 .. S
Types == {"MISS", "HOM", "HET", "UNPH"}
None == -1
NoIdx == [k \in Cl |-> None]
NaN == 3
LrVals == {1, 2, NaN}
NoUpd == [x |-> 0, before |-> [k \in Cl |-> "HOM"], idx |-> 0, thresh |-> 1]

CountUnph(t, K) == Cardinality({k \in K : t[k] = "UNPH"})
CountMiss(t, K) == Cardinality({k \in K : t[k] = "MISS"})

(* mc->unph_clusters: the unphased hets at begin_sample, in order. *)
RECURSIVE UnphFrom(_, _)
UnphFrom(t, k) == IF k >= N THEN <<>>
                  ELSE IF t[k] = "UNPH" THEN <<k>> \o UnphFrom(t, k + 1)
                  ELSE UnphFrom(t, k + 1)
UnphClusters(t) == UnphFrom(t, 0)

(* float < : false when either side is NaN. *)
FLess(a, b) == a # NaN /\ b # NaN /\ a < b

(* threshold: Arrays.sort, then the clamped rank. *)
SortedLr(l) == SortSeq(l, LAMBDA a, b : a < b)

(* Reference definitions over the environment. *)
Rank(t, k) == CountUnph(t, 0 .. k - 1)
MissRank(t, k) == CountMiss(t, 0 .. k - 1)
Flips(o) == Cardinality({i \in 1 .. Len(o) : o[i] # (IF i = 1 THEN FALSE ELSE o[i - 1])})
LastOr(o) == IF o = <<>> THEN FALSE ELSE o[Len(o)]
Max(X) == CHOOSE x \in X : \A y \in X : y <= x
(* How many times cluster k is swapped once the hets before fc are phased: the
   nearest unphased het at or before k decides, if it has been phased. *)
SwapRef(t, o, k) ==
    LET H == {j \in 0 .. k : t[j] = "UNPH"}
    IN  IF H = {} THEN 0
        ELSE LET j == Max(H)
             IN  IF j < fc /\ o[Rank(t, j) + 1] THEN 1 ELSE 0
SumOver(X, f(_)) ==
    LET sum[i \in 0 .. S] == IF i = 0 THEN 0 ELSE sum[i - 1] + (IF i \in X THEN f(i) ELSE 0)
    IN  sum[S]
(* The value at 0-based index i of the ratios in Float.compare order, by
   counting rather than sorting. *)
NthLr(l, i) == CHOOSE v \in {l[j] : j \in 1 .. Len(l)} :
                   /\ Cardinality({j \in 1 .. Len(l) : l[j] < v}) <= i
                   /\ i < Cardinality({j \in 1 .. Len(l) : l[j] <= v})
(* What update_phase leaves: the j-th unphased het stays unphased exactly
   when the j-th ratio is below the threshold. *)
KeepRef(t, l, v) ==
    LET h == UnphClusters(t)
    IN  [k \in Cl |-> IF t[k] = "UNPH"
                        THEN LET j == CHOOSE i \in 1 .. Len(h) : h[i] = k
                             IN  IF FLess(l[j], v) THEN "UNPH" ELSE "HET"
                        ELSE t[k]]

InitTypes ==
    {t \in [Cl -> Types] : \A c \in Cl : t[c] = "UNPH" => \E d \in Cl : d < c /\ t[d] = "HET"}

TypeOK ==
    /\ burnin \in BOOLEAN
    /\ types0 \in [Samples -> [Cl -> Types]] /\ types \in [Samples -> [Cl -> Types]]
    /\ snap \in [Cl -> Types]
    /\ pc \in {"begin", "bwd", "fwd", "end", "done"}
    /\ s \in 1 .. S + 1
    /\ ran \subseteq Samples
    /\ swapHaps \in BOOLEAN /\ swapCnt \in Nat
    /\ lr \in Seq(LrVals)
    /\ outs \in [Samples -> Seq(BOOLEAN)]
    /\ missIdx \in 0 .. N
    /\ bc \in -1 .. N /\ bUnph \in -1 .. N
    /\ fc \in 0 .. N /\ fUnph \in 0 .. N
    /\ hetW \in [Cl -> -1 .. N] /\ missW \in [Cl -> -1 .. N]
    /\ hetR \in [Cl -> -1 .. N] /\ missR \in [Cl -> -1 .. N]
    /\ swaps \in [Cl -> Nat]
    /\ rateSwaps \in Nat /\ rateUnph \in Nat
    /\ upd \in [x : 0 .. S, before : [Cl -> Types], idx : -1 .. N, thresh : LrVals]

Init ==
    /\ burnin \in BOOLEAN
    /\ types0 \in [Samples -> InitTypes]
    /\ types = types0
    /\ snap = types0[1]
    /\ pc = "begin" /\ s = 1
    /\ ran = {}
    /\ swapHaps = FALSE /\ swapCnt = 0
    /\ lr = <<>>
    /\ outs = [x \in Samples |-> <<>>]
    /\ missIdx = 0
    /\ bc = None /\ bUnph = None
    /\ fc = 0 /\ fUnph = 0
    /\ hetW = NoIdx /\ missW = NoIdx /\ hetR = NoIdx /\ missR = NoIdx
    /\ swaps = [k \in Cl |-> 0]
    /\ rateSwaps = 0 /\ rateUnph = 0
    /\ upd = NoUpd

(* begin_sample, and bwd_begin when it returns true. *)
Begin ==
    /\ pc = "begin" /\ s <= S
    /\ swapHaps' = FALSE /\ swapCnt' = 0
    /\ LET t == types[s]
           nUnph == CountUnph(t, Cl)
           nMiss == CountMiss(t, Cl)
       IN  IF nMiss = 0 /\ nUnph = 0
           THEN /\ pc' = "end"
                /\ UNCHANGED <<snap, ran, lr, missIdx, bc, bUnph, hetW, missW, hetR, missR, swaps>>
           ELSE /\ pc' = "bwd" /\ snap' = t /\ ran' = ran \cup {s}
                /\ lr' = <<>>
                /\ hetR' = NoIdx /\ missR' = NoIdx /\ swaps' = [k \in Cl |-> 0]
                /\ bUnph' = nUnph - 1
                /\ IF t[N - 1] = "MISS"
                     THEN missW' = [NoIdx EXCEPT ![N - 1] = nMiss - 1] /\ missIdx' = nMiss - 1
                     ELSE missW' = NoIdx /\ missIdx' = nMiss
                /\ hetW' = NoIdx
                /\ bc' = N - 2
    /\ UNCHANGED <<burnin, types0, types, s, outs, fc, fUnph, rateSwaps, rateUnph, upd>>

(* bwd_step and bwd_end_step: the values reach cluster bc, which stores them
   if it has a missing genotype, and cluster bc + 1 if it is the next
   unphased het of mc->unph_clusters. *)
BwdStep ==
    /\ pc = "bwd" /\ bc >= 0
    /\ LET hets == UnphClusters(snap)
           miss == snap[bc] = "MISS"
           het == bUnph >= 0 /\ hets[bUnph + 1] = bc + 1
       IN  /\ missW' = IF miss THEN [missW EXCEPT ![bc] = missIdx - 1] ELSE missW
           /\ missIdx' = IF miss THEN missIdx - 1 ELSE missIdx
           /\ hetW' = IF het THEN [hetW EXCEPT ![bc + 1] = bUnph] ELSE hetW
           /\ bUnph' = IF het THEN bUnph - 1 ELSE bUnph
    /\ bc' = bc - 1
    /\ UNCHANGED <<burnin, types0, types, snap, pc, s, ran, swapHaps, swapCnt, lr, outs,
                   fc, fUnph, hetR, missR, swaps, rateSwaps, rateUnph, upd>>

(* bwd_finish returns; fwd_begin. *)
FwdBegin ==
    /\ pc = "bwd" /\ bc < 0
    /\ pc' = "fwd" /\ fc' = 0 /\ fUnph' = 0
    /\ UNCHANGED <<burnin, types0, types, snap, s, ran, swapHaps, swapCnt, lr, outs, missIdx,
                   bc, bUnph, hetW, missW, hetR, missR, swaps, rateSwaps, rateUnph, upd>>

(* fwd_step at cluster fc, with phase_het and swap_haps at the next unphased
   het, then fwd_end_step with impute_alleles at a missing genotype. *)
FwdStep ==
    /\ pc = "fwd" /\ fc < N
    /\ \E swap \in BOOLEAN, v \in LrVals :
         LET c == fc
             hets == UnphClusters(snap)
             isHet == fUnph < Len(hets) /\ c = hets[fUnph + 1]
             swapEnd == IF fUnph + 1 < Len(hets) THEN hets[fUnph + 2] ELSE N
         IN  /\ isHet \/ (swap = FALSE /\ v = 1)
             /\ v = NaN => ~swap
             /\ IF isHet
                  THEN /\ hetR' = [hetR EXCEPT ![c] = fUnph]
                       /\ swapHaps' = swap
                       /\ swapCnt' = IF swap # swapHaps THEN swapCnt + 1 ELSE swapCnt
                       /\ lr' = Append(lr, v)
                       /\ outs' = [outs EXCEPT ![s] = Append(@, swap)]
                       /\ fUnph' = fUnph + 1
                       /\ swaps' = [k \in Cl |-> IF swap /\ c <= k /\ k < swapEnd
                                                 THEN swaps[k] + 1 ELSE swaps[k]]
                  ELSE UNCHANGED <<hetR, swapHaps, swapCnt, lr, outs, fUnph, swaps>>
             /\ IF snap[c] = "MISS"
                  THEN missR' = [missR EXCEPT ![c] = missIdx] /\ missIdx' = missIdx + 1
                  ELSE UNCHANGED <<missR, missIdx>>
    /\ fc' = fc + 1
    /\ UNCHANGED <<burnin, types0, types, snap, pc, s, ran, bc, bUnph, hetW, missW,
                   rateSwaps, rateUnph, upd>>

(* fwd_finish returns, or begin_sample returned false; end_sample runs
   update_phase for a sample that ran and adds its counts to the rate. *)
EndSample ==
    /\ pc = "end" \/ (pc = "fwd" /\ fc = N)
    /\ LET t == types[s]
           n == CountUnph(t, Cl)
       IN  /\ IF s \in ran /\ n > 0 /\ ~burnin
                THEN \E r \in 0 .. n, junk \in LrVals :
                       LET idx == IF r < n THEN r ELSE n - 1
                           sorted == SortedLr(lr)
                           thresh == IF 0 <= idx /\ idx < Len(sorted) THEN sorted[idx + 1] ELSE junk
                           hets == UnphClusters(t)
                       IN  /\ types' = [types EXCEPT ![s] =
                                  [k \in Cl |-> IF t[k] = "UNPH"
                                                  THEN LET j == CHOOSE i \in 1 .. Len(hets) : hets[i] = k
                                                       IN  IF FLess(lr[j], thresh) THEN "UNPH" ELSE "HET"
                                                  ELSE t[k]]]
                           /\ upd' = [x |-> s, before |-> t, idx |-> idx, thresh |-> thresh]
                ELSE UNCHANGED <<types, upd>>
           /\ rateSwaps' = rateSwaps + swapCnt
           /\ rateUnph' = rateUnph + n
    /\ pc' = "begin" /\ s' = s + 1
    /\ UNCHANGED <<burnin, types0, snap, ran, swapHaps, swapCnt, lr, outs, missIdx, bc, bUnph,
                   fc, fUnph, hetW, missW, hetR, missR, swaps>>

Finish ==
    /\ pc = "begin" /\ s > S
    /\ pc' = "done"
    /\ UNCHANGED <<burnin, types0, types, snap, s, ran, swapHaps, swapCnt, lr, outs, missIdx,
                   bc, bUnph, fc, fUnph, hetW, missW, hetR, missR, swaps, rateSwaps, rateUnph, upd>>

Stop == pc = "done" /\ UNCHANGED vars

Next == Begin \/ BwdStep \/ FwdBegin \/ FwdStep \/ EndSample \/ Finish \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(Begin \/ BwdStep \/ FwdBegin \/ FwdStep \/ EndSample \/ Finish)

-----------------------------------------------------------------------------
(* The forward cursors count the unphased hets and the missing-genotype
   clusters before c. *)
FwdLockstep ==
    pc = "fwd" => fUnph = CountUnph(snap, 0 .. fc - 1) /\ missIdx = CountMiss(snap, 0 .. fc - 1)

(* The backward cursors count down to the ones before the clusters still to
   store; the backward pass hands miss_index 0 to the forward pass. *)
BwdLockstep ==
    pc = "bwd" => /\ bUnph = CountUnph(snap, 0 .. bc + 1) - 1
                  /\ missIdx = CountMiss(snap, 0 .. bc)

(* The forward pass reads each cluster's backward values from the index the
   backward pass stored them under, its rank among clusters of its kind. *)
IndicesMatch ==
    \A k \in Cl :
        /\ hetW[k] # None => hetW[k] = Rank(snap, k)
        /\ missW[k] # None => missW[k] = MissRank(snap, k)
        /\ hetR[k] # None => hetR[k] = hetW[k] /\ hetR[k] = Rank(snap, k)
        /\ missR[k] # None => missR[k] = missW[k] /\ missR[k] = MissRank(snap, k)

(* The passes change no cluster type; only update_phase does. *)
PassKeepsTypes == pc \in {"bwd", "fwd"} => types[s] = snap

(* swap_cnt counts the changes of swap_haps, which starts each sample false. *)
SwapsCountFlips ==
    pc = "fwd" => swapCnt = Flips(outs[s]) /\ swapHaps = LastOr(outs[s])

(* Each swap covers exactly the clusters from its het to the next unphased
   het, so each cluster is swapped as the het before it decided. *)
SwapRange == pc = "fwd" => \A k \in Cl : swaps[k] = SwapRef(snap, outs[s], k)

(* A finished forward pass stored one ratio per unphased het, the ratios
   update_phase pairs with sp->unphased. *)
LrPerHet == (pc = "fwd" /\ fc = N) => Len(lr) = CountUnph(types[s], Cl)

(* The swap rate sums every phased sample's flips, and every sample's
   unphased hets before update_phase. *)
RateSums ==
    pc = "done" => /\ rateSwaps = SumOver(ran, LAMBDA x : Flips(outs[x]))
                   /\ rateUnph = SumOver(Samples, LAMBDA x : CountUnph(types0[x], Cl))

(* Nothing leaves the unphased lists in burn-in. *)
NoUpdateInBurnin == burnin => types = types0

(* threshold reads lr_sorted inside [0, n). *)
ThresholdInRange == upd.x # 0 => 0 <= upd.idx /\ upd.idx < CountUnph(upd.before, Cl)

(* update_phase keeps unphased exactly the hets whose ratio is below the
   ratio at the rank, and every other het is phased. *)
UpdateKeepsBelow ==
    (upd.x # 0 /\ pc = "begin" /\ s = upd.x + 1 /\ upd.idx \in 0 .. Len(lr) - 1) =>
        /\ upd.thresh = NthLr(lr, upd.idx)
        /\ types[upd.x] = KeepRef(upd.before, lr, upd.thresh)

(* After burn-in, every sample with an unphased het ends the run with fewer. *)
PhasingProgress ==
    (pc = "done" /\ ~burnin) =>
        \A x \in Samples : CountUnph(types0[x], Cl) > 0 => CountUnph(types[x], Cl) < CountUnph(types0[x], Cl)

(* Only UNPH ever changes, to HET. *)
Transitions ==
    [][\A x \in Samples, k \in Cl : types'[x][k] = types[x][k] \/ <<types[x][k], types'[x][k]>> = <<"UNPH", "HET">>]_vars

Terminates == <>(pc = "done")
=============================================================================
