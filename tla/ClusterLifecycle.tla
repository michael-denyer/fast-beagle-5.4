--------------------------- MODULE ClusterLifecycle ---------------------------
(***************************************************************************)
(* The clusters and marker lists of one sample across a window's runs:     *)
(* set_clusters in src/phase/sample_phase.c, which sample_phase_init runs  *)
(* once, then in each run marker_cluster_init with set_unph_clusters and   *)
(* set_has_missing (src/phase/marker_cluster.c), and update_phase          *)
(* (src/phase/phase_baum1.c) shrinking sp->unphased after burn-in. Markers *)
(* are 0 .. M-1. Init chooses each marker's genotype (homozygous, het, or  *)
(* missing) and the stage-1 overlap start; sp->missing lists the missing   *)
(* markers and sp->unphased the hets at or after the overlap start other   *)
(* than the sample's first het, as PbwtPhaser.indices (pbwt_phaser.c)      *)
(* builds them. MaxSize stands for the 255-marker cluster limit. At each   *)
(* marker the environment chooses whether it lies beyond 0.005 cM from the *)
(* cluster start. Init chooses how many of the Runs runs are burn-in, at   *)
(* least one of each. update_phase keeps the j-th unphased marker when     *)
(* the j-th likelihood ratio is below the threshold; which ones it keeps   *)
(* is the environment's choice (HmmPass checks the threshold). Each step   *)
(* is one pass of a set_clusters, set_unph_clusters or set_has_missing     *)
(* loop, or one update_phase.                                              *)
(***************************************************************************)
EXTENDS Integers, Sequences, FiniteSets

CONSTANTS M, MaxSize, Runs

ASSUME M \in Nat \ {0} /\ MaxSize \in Nat \ {0} /\ Runs \in Nat \ {0, 1}

VARIABLES
    kind,      \* per marker: "HOM", "HET" or "MISS"
    unph0,     \* sp->unphased at sample_phase_init
    unph,      \* sp->unphased
    prevUnph,  \* sp->unphased when the current run started
    nBurn,     \* the burn-in runs
    pc,        \* "clust", "unph", "miss", "update" or "done"
    run,       \* runs started
    m, lastEnd, prevMH,  \* set_clusters: m, last_end, prev_is_missing_or_het
    sizes,     \* sp->clust_size
    j,         \* the cluster loop's index in set_unph_clusters and set_has_missing
    ui, nextU, \* set_unph_clusters: unph_index and next_unph
    unphCl,    \* mc->unph_clusters
    mi,        \* set_has_missing: i
    hasMiss    \* mc->has_missing

vars == <<kind, unph0, unph, prevUnph, nBurn, pc, run, m, lastEnd, prevMH, sizes,
          j, ui, nextU, unphCl, mi, hasMiss>>

Mk == 0 .. M - 1
Kinds == {"HOM", "HET", "MISS"}
Inf == M + 1
Range(q) == {q[i] : i \in 1 .. Len(q)}

RECURSIVE ListFrom(_, _)
ListFrom(S, k) == IF k >= M THEN <<>> ELSE IF k \in S THEN <<k>> \o ListFrom(S, k + 1) ELSE ListFrom(S, k + 1)
List(S) == ListFrom(S, 0)

Missing == List({k \in Mk : kind[k] = "MISS"})
Min(X) == CHOOSE x \in X : \A y \in X : x <= y

(* PbwtPhaser.indices: the hets at or after ov, other than the first het. *)
Indices(kd, ov) ==
    LET hets == {k \in Mk : kd[k] = "HET"}
    IN  List({k \in hets : k >= ov /\ k # Min(hets)})

Sum(q, n) == LET f[i \in 0 .. n] == IF i = 0 THEN 0 ELSE f[i - 1] + q[i] IN f[n]
NClust == Len(sizes)
Ends == [c \in 0 .. NClust - 1 |-> Sum(sizes, c + 1)]
Start(c) == IF c = 0 THEN 0 ELSE Ends[c - 1]
ClusterOf(k) == CHOOSE c \in 0 .. NClust - 1 : Start(c) <= k /\ k < Ends[c]

(* unph_index < n_unph ? sp->unphased[unph_index++] : INT_MAX, as
   <<next_unph, unph_index>>. *)
FetchIn(u, i) == IF i < Len(u) THEN <<u[i + 1], i + 1>> ELSE <<Inf, i>>
Fetch(i) == FetchIn(unph, i)
(* The while loop: next_unph past the markers before end. *)
RECURSIVE While(_, _)
While(p, end) == IF p[1] < end THEN While(Fetch(p[2]), end) ELSE p
(* The for loop over sp->missing below end. *)
RECURSIVE MissSkip(_, _)
MissSkip(i, end) == IF i < Len(Missing) /\ Missing[i + 1] < end THEN MissSkip(i + 1, end) ELSE i

IsSubSeq(a, b) == \E f \in [1 .. Len(a) -> 1 .. Len(b)] :
                      /\ \A i \in 1 .. Len(a) : a[i] = b[f[i]]
                      /\ \A i \in 1 .. Len(a) - 1 : f[i] < f[i + 1]

TypeOK ==
    /\ kind \in [Mk -> Kinds]
    /\ unph \in Seq(Mk) /\ unph0 \in Seq(Mk) /\ prevUnph \in Seq(Mk)
    /\ nBurn \in 1 .. Runs - 1
    /\ pc \in {"clust", "unph", "miss", "update", "done"}
    /\ run \in 0 .. Runs
    /\ m \in 0 .. M /\ lastEnd \in 0 .. M /\ prevMH \in BOOLEAN
    /\ sizes \in Seq(1 .. MaxSize)
    /\ j \in 0 .. M /\ ui \in 0 .. M /\ nextU \in 0 .. Inf
    /\ unphCl \in Seq(0 .. M - 1)
    /\ mi \in 0 .. M
    /\ hasMiss \in Seq(BOOLEAN)

Init ==
    /\ kind \in [Mk -> Kinds]
    /\ \E ov \in 0 .. M : unph0 = Indices(kind, ov)
    /\ unph = unph0 /\ prevUnph = unph0
    /\ nBurn \in 1 .. Runs - 1
    /\ pc = "clust" /\ run = 0
    /\ m = 0 /\ lastEnd = 0 /\ prevMH = FALSE /\ sizes = <<>>
    /\ j = 0 /\ ui = 0 /\ nextU = Inf /\ unphCl = <<>> /\ mi = 0 /\ hasMiss = <<>>

(* One pass of the set_clusters loop at marker m. *)
ClustStep ==
    /\ pc = "clust" /\ m < M
    /\ \E far \in BOOLEAN :
         LET size == m - lastEnd
             mh == kind[m] = "MISS" \/ kind[m] = "HET"
             cut == prevMH \/ mh \/ far \/ size = MaxSize
         IN  /\ IF cut /\ m > 0
                  THEN sizes' = Append(sizes, size) /\ lastEnd' = m
                  ELSE UNCHANGED <<sizes, lastEnd>>
             /\ prevMH' = mh
    /\ m' = m + 1
    /\ UNCHANGED <<kind, unph0, unph, prevUnph, nBurn, pc, run, j, ui, nextU, unphCl, mi, hasMiss>>

(* marker_cluster_init up to its loops. *)
BeginRun(u) ==
    /\ pc' = "unph" /\ run' = run + 1 /\ prevUnph' = u
    /\ j' = 0 /\ unphCl' = <<>> /\ mi' = 0 /\ hasMiss' = <<>>
    /\ nextU' = FetchIn(u, 0)[1] /\ ui' = FetchIn(u, 0)[2]

(* The last cluster, then the first run. *)
ClustEnd ==
    /\ pc = "clust" /\ m = M
    /\ sizes' = Append(sizes, M - lastEnd)
    /\ BeginRun(unph)
    /\ UNCHANGED <<kind, unph0, unph, nBurn, m, lastEnd, prevMH>>

(* One pass of the set_unph_clusters loop at cluster j. *)
UnphStep ==
    /\ pc = "unph"
    /\ IF j < NClust
         THEN LET end == Ends[j]
              IN  /\ IF nextU < end
                       THEN LET q == While(Fetch(ui), end)
                            IN  unphCl' = Append(unphCl, j) /\ nextU' = q[1] /\ ui' = q[2]
                       ELSE UNCHANGED <<unphCl, nextU, ui>>
                  /\ j' = j + 1 /\ UNCHANGED pc
         ELSE pc' = "miss" /\ j' = 0 /\ UNCHANGED <<unphCl, nextU, ui>>
    /\ UNCHANGED <<kind, unph0, unph, prevUnph, nBurn, run, m, lastEnd, prevMH, sizes, mi, hasMiss>>

(* One pass of the set_has_missing loop at cluster j. *)
MissStep ==
    /\ pc = "miss"
    /\ IF j < NClust
         THEN LET i2 == MissSkip(mi, Ends[j])
              IN  hasMiss' = Append(hasMiss, i2 > mi) /\ mi' = i2 /\ j' = j + 1 /\ UNCHANGED pc
         ELSE pc' = "update" /\ UNCHANGED <<hasMiss, mi, j>>
    /\ UNCHANGED <<kind, unph0, unph, prevUnph, nBurn, run, m, lastEnd, prevMH, sizes, ui, nextU, unphCl>>

(* end_sample's update_phase after the passes: after burn-in it keeps the
   unphased markers whose ratio passed, in order. Then the next run, or the
   window's end. *)
Update ==
    /\ pc = "update"
    /\ IF Len(unph) > 0 /\ run > nBurn
         THEN \E keep \in SUBSET (1 .. Len(unph)) :
                unph' = [i \in 1 .. Cardinality(keep) |->
                            unph[CHOOSE x \in keep : Cardinality({y \in keep : y < x}) = i - 1]]
         ELSE UNCHANGED unph
    /\ IF run < Runs
         THEN BeginRun(unph')
         ELSE pc' = "done" /\ UNCHANGED <<run, prevUnph, j, ui, nextU, unphCl, mi, hasMiss>>
    /\ UNCHANGED <<kind, unph0, nBurn, m, lastEnd, prevMH, sizes>>

Stop == pc = "done" /\ UNCHANGED vars

Next == ClustStep \/ ClustEnd \/ UnphStep \/ MissStep \/ Update \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(ClustStep \/ ClustEnd \/ UnphStep \/ MissStep \/ Update)

-----------------------------------------------------------------------------
(* The clusters partition the markers, and every het or missing genotype is
   a cluster of its own; TypeOK bounds each size by MaxSize. *)
Partition ==
    pc # "clust" =>
        /\ Sum(sizes, NClust) = M
        /\ \A k \in Mk : kind[k] # "HOM" => sizes[ClusterOf(k) + 1] = 1

(* mc->unph_clusters lists, in order, the clusters that hold an unphased
   marker. *)
UnphClustersMatch ==
    pc \in {"miss", "update"} => unphCl = List({ClusterOf(k) : k \in Range(unph)})

(* The j-th unphased het has the j-th cluster, so phase_het fills one ratio
   per entry update_phase reads, each for its own marker. *)
LrPairsWithMarker ==
    pc \in {"miss", "update"} =>
        /\ Len(unphCl) = Len(unph)
        /\ \A i \in 1 .. Len(unph) : unphCl[i] = ClusterOf(unph[i])

(* No unphased het is in cluster 0, the one bwd_end_step never stores. *)
NoUnphAtFirst == pc \in {"miss", "update"} => \A i \in 1 .. Len(unphCl) : unphCl[i] > 0

(* mc->has_missing marks the clusters of the missing markers, one each. *)
HasMissingMatch ==
    pc = "update" =>
        /\ Len(hasMiss) = NClust
        /\ \A c \in 0 .. NClust - 1 : hasMiss[c + 1] <=> \E k \in Range(Missing) : ClusterOf(k) = c
        /\ Cardinality({c \in 1 .. NClust : hasMiss[c]}) = Len(Missing)

(* sp->unphased stays increasing, and each update keeps a sublist of it. *)
UnphasedShrinks ==
    /\ \A i \in 1 .. Len(unph) - 1 : unph[i] < unph[i + 1]
    /\ IsSubSeq(unph, prevUnph) /\ IsSubSeq(unph, unph0)

(* No update_phase shrinks sp->unphased in a burn-in run. *)
BurninKeepsList == [][pc = "update" /\ run <= nBurn => unph' = unph]_vars

Terminates == <>(pc = "done")
=============================================================================
