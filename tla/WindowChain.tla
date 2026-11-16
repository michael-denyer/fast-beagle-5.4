------------------------------ MODULE WindowChain ------------------------------
(***************************************************************************)
(* The chain of windows the reader thread builds in                        *)
(* src/vcf/sliding_window.c: read_windows, next_end_cm, next_targ_window   *)
(* (Java's TargSlidingWindow.Reader), next_ref_window                      *)
(* (RefTargSlidingWindow.Reader) and overlap_start. The hand-over of each  *)
(* window to the caller is tla/SlidingWindow.tla; here a window is one     *)
(* step. Mode "targ" has no reference; "ref" and "impute" read the         *)
(* reference with impute=false and impute=true. Init chooses the target    *)
(* records T and reference records R, each a marker [c (chromosome), p     *)
(* (position), a (allele set)], sorted by chromosome and position. In      *)
(* "targ" mode records at one position stand for markers with different    *)
(* alleles, so a is 1 and only c and p are used. The genetic map is the    *)
(* identity, one cM per base, so window=W and overlap=O are in bases; this *)
(* hides the rounding of genetic_map_base_pos and flat stretches of a real *)
(* map. As in Beagle 5.4, a window that starts with an overlap ends W - O  *)
(* past the end of the window before, so it can read no new record, and    *)
(* there is no window-markers cap, no overlap-markers bound on the overlap *)
(* search and no check that a window advances. empty_window_error and      *)
(* Java's empty-window exception end the chain at pc = "empty".            *)
(***************************************************************************)
EXTENDS Integers, Sequences, FiniteSets

CONSTANTS Mode, NT, NR, NC, P, A, W, O

ASSUME /\ Mode \in {"targ", "ref", "impute"}
       /\ NT \in Nat /\ NR \in Nat /\ (Mode = "targ" => NR = 0 /\ A = 1)
       /\ {NC, P, A} \subseteq Nat \ {0}
       /\ O \in Nat /\ W > O /\ 10 * W >= 11 * O

Markers == [c : 1..NC, p : 1..P, a : 1..A]

Before(m, n) == m.c < n.c \/ (m.c = n.c /\ m.p <= n.p)

(* Input files sorted by chromosome and position; a reference or target VCF
   holds a marker once. *)
Inputs(n) == {s \in [1..n -> Markers] :
                /\ \A j \in 1..n - 1 : Before(s[j], s[j + 1])
                /\ Mode # "targ" => \A j, k \in 1..n : j # k => s[j] # s[k]}

Entries == [r : 1..NR, in : BOOLEAN]

(* A bound on how many windows in a row start with an overlap: each reads a
   record or moves the end past one more position. *)
MaxRun == NT + NR + P

VARIABLES
    T, R,     \* the target and reference records
    pc,       \* "start", "stop", "empty" or "nodata"
    next,     \* the index of sw->next_targ; NT + 1 for NULL
    rnext,    \* the index of sw->next_ref; NR + 1 for NULL
    endCm,    \* sw->end_cm
    wc,       \* the chromosome the last window read
    tov,      \* sw->targ_overlap, as indices into T
    rov,      \* sw->ref_overlap, as indices into R with in_targ
    twin,     \* the target records of the last window made
    rwin,     \* the reference records of the last window made
    inLen,    \* how many records of the last window were carried into it
    seen,     \* the target records of every window made
    done,     \* sw->done
    anchor,   \* the position next_end_cm last started from with an empty overlap
    run       \* how many windows since then started with an overlap

vars == <<T, R, pc, next, rnext, endCm, wc, tov, rov, twin, rwin, inLen, seen, done, anchor, run>>

Range(s) == {s[j] : j \in 1..Len(s)}

TypeOK ==
    /\ T \in Inputs(NT) /\ R \in Inputs(NR)
    /\ pc \in {"start", "stop", "empty", "nodata"}
    /\ next \in 1..NT + 1 /\ rnext \in 1..NR + 1
    /\ endCm \in 0..P + W + MaxRun * (W - O) /\ wc \in 0..NC
    /\ tov \in Seq(1..NT) /\ rov \in Seq(Entries)
    /\ twin \in Seq(1..NT) /\ rwin \in Seq(Entries)
    /\ inLen \in Nat /\ seen \subseteq 1..NT
    /\ done \in BOOLEAN
    /\ anchor \in 0..P /\ run \in 0..MaxRun

(* read_windows reads the first record of each input. *)
Init ==
    /\ T \in Inputs(NT) /\ R \in Inputs(NR)
    /\ next = 1 /\ rnext = 1
    /\ endCm = 0 /\ wc = 0
    /\ pc = IF NT = 0 \/ (Mode # "targ" /\ NR = 0) THEN "nodata" ELSE "start"
    /\ tov = <<>> /\ rov = <<>> /\ twin = <<>> /\ rwin = <<>>
    /\ inLen = 0 /\ seen = {} /\ done = FALSE /\ anchor = 0 /\ run = 0

(* next_end_cm: window= past the record at position p when the overlap is
   empty, else window-overlap past the previous end. *)
NextEnd(empty, p) == IF empty THEN p + W ELSE endCm + (W - O)

(* run stops at MaxRun + 1, past the bound TypeOK holds it to, so a chain
   that repeats a window forever repeats a state. *)
Restart(empty, p) == /\ anchor' = IF empty THEN p ELSE anchor
                     /\ run' = IF empty THEN 0 ELSE IF run > MaxRun THEN run ELSE run + 1

(* overlap_start over the positions ps of a window's markers: binary search
   for the position overlap= before endGen, then back to the first marker at
   that position. Indices are 0-based, as in the C. *)
RECURSIVE FirstWithPos(_, _)
FirstWithPos(ps, j) == IF j > 0 /\ ps[j] = ps[j + 1] THEN FirstWithPos(ps, j - 1) ELSE j

RECURSIVE Search(_, _, _, _)
Search(ps, key, low, high) ==
    IF low > high THEN FirstWithPos(ps, IF high > 0 THEN high ELSE 0)
    ELSE LET mid == (low + high) \div 2
         IN  IF ps[mid + 1] < key THEN Search(ps, key, mid + 1, high)
             ELSE IF ps[mid + 1] > key THEN Search(ps, key, low, mid - 1)
             ELSE FirstWithPos(ps, mid)

OverlapStart(ps, endGen) == Search(ps, endGen - O, 0, Len(ps) - 1)

(* The read loop of next_targ_window. *)
RECURSIVE TLoop(_, _, _, _)
TLoop(recs, k, chrom, endPos) ==
    IF k <= NT /\ T[k].c = chrom /\ T[k].p < endPos
    THEN TLoop(Append(recs, k), k + 1, chrom, endPos)
    ELSE <<recs, k>>

(* next_targ_window and the carried overlap. The stages are bound with
   \E x \in {e}, not LET, because TLC evaluates a LET again each time it
   checks the fairness condition; with LET the reference modes ran ten times
   slower. *)
TargWindow ==
    /\ pc = "start" /\ ~done /\ Mode = "targ"
    /\ endCm' = NextEnd(tov = <<>>, T[next].p)
    /\ Restart(tov = <<>>, T[next].p)
    /\ wc' = T[next].c
    /\ LET chrom == T[next].c
       IN  \E lp \in {TLoop(tov, next, chrom, endCm')} :
           LET recs == lp[1]
               k == lp[2]
               last == k > NT
               chromEnd == last \/ T[k].c # chrom
           IN  IF recs = <<>>
               THEN pc' = "empty" /\ UNCHANGED <<next, tov, twin, inLen, seen, done>>
               ELSE \E ovStart \in {IF chromEnd THEN Len(recs)
                                     ELSE OverlapStart([j \in 1..Len(recs) |-> T[recs[j]].p], T[recs[Len(recs)]].p)} :
                    /\ next' = k
                    /\ tov' = SubSeq(recs, ovStart + 1, Len(recs))
                    /\ twin' = recs /\ inLen' = Len(tov) /\ seen' = seen \cup Range(recs)
                    /\ done' = last
                    /\ UNCHANGED pc
    /\ UNCHANGED <<T, R, rnext, rov, rwin>>

Impute == Mode = "impute"

(* advanceRefItToChrom: it stops at the last reference record. *)
RECURSIVE Advance(_, _)
Advance(rk, chrom) == IF R[rk].c # chrom /\ rk < NR THEN Advance(rk + 1, chrom) ELSE rk

(* The inner loop of next_ref_window over the reference records before
   target record k. *)
RECURSIVE Skip(_, _, _, _)
Skip(ref, rk, k, chrom) ==
    IF rk <= NR /\ R[rk].c = chrom /\ (R[rk].p < T[k].p \/ (R[rk].p = T[k].p /\ R[rk] # T[k]))
    THEN Skip(IF Impute THEN Append(ref, [r |-> rk, in |-> FALSE]) ELSE ref, rk + 1, k, chrom)
    ELSE [ref |-> ref, rk |-> rk]

(* The target loop of next_ref_window; s is the inner loop's result for
   target record k. *)
RECURSIVE RLoop(_, _, _, _, _, _), RMatch(_, _, _, _, _)
RLoop(targ, ref, k, rk, chrom, endPos) ==
    IF k <= NT /\ T[k].c = chrom /\ T[k].p < endPos
    THEN RMatch(Skip(ref, rk, k, chrom), targ, k, chrom, endPos)
    ELSE [targ |-> targ, ref |-> ref, k |-> k, rk |-> rk]

RMatch(s, targ, k, chrom, endPos) ==
    IF s.rk <= NR /\ R[s.rk] = T[k]
    THEN RLoop(Append(targ, k), Append(s.ref, [r |-> s.rk, in |-> TRUE]), k + 1, s.rk + 1, chrom, endPos)
    ELSE RLoop(targ, s.ref, k + 1, s.rk, chrom, endPos)

(* The impute=true loop over the reference records after the last target. *)
RECURSIVE Trail(_, _, _, _)
Trail(ref, rk, chrom, endPos) ==
    IF Impute /\ rk <= NR /\ R[rk].c = chrom /\ R[rk].p < endPos
    THEN Trail(Append(ref, [r |-> rk, in |-> FALSE]), rk + 1, chrom, endPos)
    ELSE [ref |-> ref, rk |-> rk]

(* next_ref_window and the carried overlaps. next_end_cm reads the reference
   record advanceRefItToChrom stopped at. MarkerIndices.targOverlapStart is
   the first target marker at or after the reference overlap start. *)
RefWindow ==
    /\ pc = "start" /\ ~done /\ Mode # "targ"
    /\ wc' = T[next].c
    /\ LET chrom == T[next].c
       IN  \E rk0 \in {Advance(rnext, chrom)} :
           /\ endCm' = NextEnd(rov = <<>>, R[rk0].p)
           /\ Restart(rov = <<>>, R[rk0].p)
           /\ \E lp \in {RLoop(tov, rov, next, rk0, chrom, endCm')} :
              \E tr \in {Trail(lp.ref, lp.rk, chrom, endCm')} :
              LET targ == lp.targ
                  ref == tr.ref
                  k == lp.k
                  rk == tr.rk
                  last == k > NT \/ rk > NR
                  chromEnd == last \/ R[ref[1].r].c # R[rk].c
              IN  IF targ = <<>> \/ ref = <<>>
                  THEN pc' = "empty" /\ UNCHANGED <<next, rnext, tov, rov, twin, rwin, inLen, done>>
                  ELSE \E ovStart \in {IF chromEnd THEN Len(ref)
                                        ELSE OverlapStart([j \in 1..Len(ref) |-> R[ref[j].r].p], endCm' - 1)} :
                       LET tStart == Cardinality({j \in 1..ovStart : ref[j].in})
                       IN  /\ next' = k /\ rnext' = rk
                           /\ tov' = SubSeq(targ, tStart + 1, Len(targ))
                           /\ rov' = SubSeq(ref, ovStart + 1, Len(ref))
                           /\ twin' = targ /\ rwin' = ref /\ inLen' = Len(rov)
                           /\ done' = last
                           /\ UNCHANGED pc
    /\ UNCHANGED <<T, R, seen>>

(* read_windows hands over NULL after the last window. *)
Finish == pc = "start" /\ done /\ pc' = "stop"
          /\ UNCHANGED <<T, R, next, rnext, endCm, wc, tov, rov, twin, rwin, inLen, seen, done, anchor, run>>

Stop == pc # "start" /\ UNCHANGED vars

Next == TargWindow \/ RefWindow \/ Finish \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(TargWindow \/ RefWindow \/ Finish)

-----------------------------------------------------------------------------
Made == twin # <<>>

IsSuffix(s, t) == Len(s) <= Len(t) /\ s = SubSeq(t, Len(t) - Len(s) + 1, Len(t))

(* A window is a run of consecutive target records. *)
WinInterval == \A j \in 1..Len(twin) - 1 : twin[j + 1] = twin[j] + 1

(* The windows made so far cover a prefix of the target records, and all of
   them once the chain ends. *)
Covered == /\ \E m \in 0..NT : seen = 1..m
           /\ pc = "stop" => seen = 1..NT

(* A window lies on one chromosome. *)
OneChrom == \A a, b \in Range(twin) : T[a].c = T[b].c

(* The overlap carried into the next window is a suffix of the window made. *)
CarrySuffix == pc = "start" /\ Made => IsSuffix(tov, twin) /\ IsSuffix(rov, rwin)

(* Nothing is carried across a chromosome boundary or past the end: of the
   target in "targ" mode, of the reference otherwise. *)
BoundaryEmpty ==
    pc = "start" /\ Made =>
    IF Mode = "targ" THEN next > NT \/ T[next].c # T[twin[1]].c => tov = <<>>
    ELSE done \/ R[rnext].c # R[rwin[1].r].c => tov = <<>> /\ rov = <<>>

(* done is set by the window holding the last target record and no other. *)
DoneAtLast == Made => (done <=> NT \in Range(twin))

IsIn(e) == e.in

(* The in-target reference records pair, in order, with target records of
   the same marker: in the window made and in what it carries. *)
Pairs(rs, ts) == LET m == SelectSeq(rs, IsIn)
                 IN  Len(m) = Len(ts) /\ \A j \in 1..Len(ts) : R[m[j].r] = T[ts[j]]

Paired == Made => Pairs(rwin, twin) /\ Pairs(rov, tov)

(* The reference records a window read lie on one chromosome, and so do the
   ones carried into it: a window spans two chromosomes only when it starts
   with an overlap from the chromosome before (Window.chromIndex). *)
ChromsBesideCarry ==
    Mode # "targ" =>
    /\ \A a, b \in inLen + 1..Len(rwin) : R[rwin[a].r].c = R[rwin[b].r].c
    /\ \A a, b \in 1..inLen : R[rwin[a].r].c = R[rwin[b].r].c

(* A window holds every record of its chromosome before its end: the read
   loops stop at the end, not earlier. *)
ReadsToEnd ==
    pc = "start" /\ Made =>
    /\ next > NT \/ T[next].c # wc \/ T[next].p >= endCm
    /\ Impute => rnext > NR \/ R[rnext].c # wc \/ R[rnext].p >= endCm

Max(S) == CHOOSE x \in S : \A y \in S : y <= x

(* How many of the sorted positions ps lie at or after the last position at
   or before key, or all of them when none is: the overlap overlap_start
   finds, searching the whole window. *)
CarryLen(ps, key) == LET below == {j \in 1..Len(ps) : ps[j] <= key}
                         q == IF below = {} THEN ps[1] ELSE ps[Max(below)]
                     IN  Cardinality({j \in 1..Len(ps) : ps[j] >= q})

(* The overlap carried to a window on the same chromosome starts overlap=
   before the window's last target marker (TargSlidingWindow) or before
   its end (RefTargSlidingWindow). *)
OverlapFromEnd ==
    pc = "start" /\ Made =>
    IF Mode = "targ"
    THEN next <= NT /\ T[next].c = wc =>
         Len(tov) = CarryLen([j \in 1..Len(twin) |-> T[twin[j]].p], T[twin[Len(twin)]].p - O)
    ELSE ~done /\ R[rwin[1].r].c = R[rnext].c =>
         Len(rov) = CarryLen([j \in 1..Len(rwin) |-> R[rwin[j].r].p], endCm - 1 - O)

(* Beagle 5.4's windows end at fixed steps: the k-th window after one that
   started without an overlap ends window + k * (window - overlap) past the
   record that window's end was measured from. *)
EndOnGrid == pc = "start" /\ Made => endCm = anchor + W + run * (W - O)

(* The chain ends: the last window or an error. *)
Terminates == <>(pc # "start")
=============================================================================
