----------------------------- MODULE RefSeqGroup -----------------------------
(***************************************************************************)
(* The sequence-coding buffer of src/vcf/ref_it.c (fill at line 133,       *)
(* add_rec 103, apply_seq_coding 92, flush_compressed 70, ref_it_next 179  *)
(* and the fill in ref_it_open 169) with seq_coder3_add and                *)
(* seq_coder3_flush of src/bref/seq_coder3.c. fill codes whole blocks of B *)
(* records from block_reader_next, the first one record longer, until a    *)
(* flush has made records ready, and flushes when a block comes back       *)
(* empty, as block_reader_next returns NULL after the last record and on   *)
(* every later call. A record that passes the marker filter goes to        *)
(* low_freq, or to the coder with a NULL placeholder in low_freq. A new    *)
(* chromosome, a full coder (add returns false, then the add is retried    *)
(* once), and the end of the input flush the coder: its records fill the   *)
(* placeholders in order and low_freq moves to ready. The coder refuses a  *)
(* record once its group would carry MaxNSeq allele strings                *)
(* (seq_coder3.c:131, as SeqCoder3.setAlleleMap tests size() >= maxNSeq),  *)
(* so a record that carries MaxNSeq strings by itself is refused by the    *)
(* retry too. Beagle 5.4's RefIt only asserts that the retry succeeds, so  *)
(* without -ea the record is dropped, its placeholder stays, and the next  *)
(* flush throws IndexOutOfBoundsException (ref_it.c:79-82).                *)
(*                                                                         *)
(* A record is a function from the H haplotypes to alleles, written with   *)
(* the alleles numbered in order of first appearance; the coder depends    *)
(* only on which haplotypes share an allele. Init chooses each record's    *)
(* chromosome, whether the filter drops it, and whether it is sequence     *)
(* coded. apply_seq_coding's bound on n_alleles is kept, so a coded record *)
(* has at most MaxCoded alleles; its bound on the major allele's count is  *)
(* a free choice, which gives the coder more inputs than the code does.    *)
(* The major allele is any allele of greatest count. The coder's sequence  *)
(* indices are modeled up to renaming: hapSeq names each sequence by its   *)
(* smallest haplotype. Each step is one record read, one record filtered   *)
(* with its flush decision, one add, one block boundary, or one call.      *)
(***************************************************************************)
EXTENDS Naturals, Sequences, FiniteSets

CONSTANTS N, H, MaxNSeq, C, B

ASSUME MaxNSeq >= 2 /\ B >= 1

Haps == 1..H

(* it->max_seq_coded_alleles in ref_it_open. *)
MaxCoded == IF MaxNSeq < 255 THEN MaxNSeq ELSE 255

Range(f) == {f[x] : x \in DOMAIN f}
Count(v, a) == Cardinality({h \in Haps : v[h] = a})

(* Allele vectors with alleles numbered in order of first appearance. *)
Vectors == {v \in [Haps -> 0..(MaxCoded - 1)] :
                \A h \in Haps : v[h] > 0 => \E f \in 1..(h - 1) : v[f] = v[h] - 1}
MostFrequent(v) == {a \in Range(v) : \A b \in Range(v) : Count(v, b) <= Count(v, a)}

Zero == [h \in Haps |-> 0]
Records ==
    {[kind |-> "excl", chrom |-> 1, v |-> Zero, maj |-> 0]}
    \cup {[kind |-> "low", chrom |-> c, v |-> Zero, maj |-> 0] : c \in 1..C}
    \cup UNION {{[kind |-> "coded", chrom |-> c, v |-> v, maj |-> m] : m \in MostFrequent(v)} :
                    c \in 1..C, v \in Vectors}

Fresh == [h \in Haps |-> 1]

Item(id, g) == [k |-> "rec", id |-> id, grp |-> g]
End == [k |-> "end", id |-> 0, grp |-> {}]
Throw == [k |-> "throw", id |-> 0, grp |-> {}]
None == [k |-> "none", id |-> 0, grp |-> {}]

VARIABLES
    Recs,      \* the records block_reader_next returns, in file order
    i,         \* the record block_reader_next returns next; N + 1 for NULL
    nblk,      \* it->n_blocks
    blk,       \* recs[0..n) in fill, as record numbers
    bj,        \* the records of blk that add_rec has taken
    cur,       \* the record add_rec is placing
    lastChrom, \* it->last_chrom; 0 for -1, none yet
    lowFreq,   \* it->low_freq: record numbers, 0 for a NULL placeholder
    group,     \* the coder's records, c->recs
    hapSeq,    \* c->hap_to_seq, each sequence named by its smallest haplotype
    nseq,      \* c->n_seq_map
    ready,     \* it->ready
    pending,   \* the record ref_it_next returns once fill returns
    dropped,   \* a refused retry left a placeholder with no record
    pc,        \* "fill", "read", "add", "code", "retry", "ret", "idle" or "throw"
    out,       \* what ref_it_next returned, in order
    expected   \* Expected for these records, fixed by Init

vars == <<Recs, i, nblk, blk, bj, cur, lastChrom, lowFreq, group, hapSeq, nseq, ready,
          pending, dropped, pc, out, expected>>

(* The number of distinct allele strings the haplotypes carry over the
   records of g. *)
NSeq(g) == Cardinality({[k \in 1..Len(g) |-> Recs[g[k]].v[h]] : h \in Haps})

-----------------------------------------------------------------------------
(* The reference: Java 5.4 RefIt over these records. The flushes, in order:
   each has the position of the record it comes before (N + 1 for the
   end), the coded group it hands out, and whether it throws. A coded
   record joins the open group unless the chromosome changed or the group
   would then carry MaxNSeq allele strings; a record that carries MaxNSeq
   by itself is dropped, and the next flush throws. *)
RECURSIVE Flushes(_, _, _, _)
Flushes(j, g, c, d) ==
    LET F(p) == [pos |-> p, grp |-> g, throws |-> d]
    IN  IF j > N THEN <<F(N + 1)>>
        ELSE LET r == Recs[j] IN
             IF r.kind = "excl" THEN Flushes(j + 1, g, c, d)
             ELSE IF c # 0 /\ r.chrom # c
               THEN IF d THEN <<F(j)>> ELSE <<F(j)>> \o Flushes(j, <<>>, r.chrom, d)
             ELSE IF r.kind = "low" THEN Flushes(j + 1, g, r.chrom, d)
             ELSE IF NSeq(Append(g, j)) < MaxNSeq THEN Flushes(j + 1, Append(g, j), r.chrom, d)
             ELSE IF d THEN <<F(j)>>
             ELSE <<F(j)>> \o
                  IF NSeq(<<j>>) < MaxNSeq THEN Flushes(j + 1, <<j>>, r.chrom, FALSE)
                  ELSE Flushes(j + 1, <<>>, r.chrom, TRUE)

(* BlockLineReader's blocks: RefIt.combine puts the first data line in
   front of the first block. The end's flush comes in the empty block
   after the last. *)
Blk(p) == IF p <= B + 1 THEN 1 ELSE (p - B - 2) \div B + 2
FlushBlk(p) == IF p <= N THEN Blk(p) ELSE IF N = 0 THEN 1 ELSE Blk(N) + 1

(* Every kept record in file order, each coded record with its group, and
   then NULL. When a flush throws, the caller has seen the records that
   flushes in earlier blocks made ready, since fillRecBuffer codes whole
   blocks before it returns, except the last: next() had taken it and
   called fillRecBuffer, which threw. *)
Expected ==
    LET ev == Flushes(1, <<>>, 0, FALSE)
        last == ev[Len(ev)]
        EvOf(p) == CHOOSE k \in 1..Len(ev) : ev[k].pos > p /\ \A m \in 1..(k - 1) : ev[m].pos <= p
        Seen(p) == /\ Recs[p].kind # "excl"
                   /\ p < last.pos
                   /\ last.throws => FlushBlk(ev[EvOf(p)].pos) < FlushBlk(last.pos)
        shown == SelectSeq([p \in 1..N |-> p], Seen)
        ItemOf(p) == Item(p, IF Recs[p].kind = "coded" THEN Range(ev[EvOf(p)].grp) ELSE {})
    IN  IF last.throws
        THEN [m \in 1..(Len(shown) - 1) |-> ItemOf(shown[m])] \o <<Throw>>
        ELSE [m \in 1..Len(shown) |-> ItemOf(shown[m])] \o <<End>>

-----------------------------------------------------------------------------
Seqs == Range(hapSeq)
Placeholders == Cardinality({p \in 1..Len(lowFreq) : lowFreq[p] = 0})

(* set_allele_map for record r: the new sequences its non-major alleles
   push, after the first per sequence, and add_major_allele's push for a
   sequence that holds both the major allele and another. *)
NewNSeq(r) ==
    LET NonMajor(s) == {r.v[h] : h \in {g \in Haps : hapSeq[g] = s /\ r.v[g] # r.maj}}
        HasMajor(s) == \E h \in Haps : hapSeq[h] = s /\ r.v[h] = r.maj
    IN  nseq
        + Cardinality(UNION {{<<s, a>> : a \in NonMajor(s)} : s \in Seqs})
        - Cardinality({s \in Seqs : NonMajor(s) # {}})
        + Cardinality({s \in Seqs : HasMajor(s) /\ NonMajor(s) # {}})

(* seq_coder3_add's new hap_to_seq: haplotypes stay together while they
   share the old sequence and the record's allele. *)
Refine(v) ==
    [h \in Haps |-> CHOOSE g \in Haps :
        /\ hapSeq[g] = hapSeq[h] /\ v[g] = v[h]
        /\ \A f \in Haps : (hapSeq[f] = hapSeq[h] /\ v[f] = v[h]) => g <= f]

(* seq_coder3_add succeeds for record cur and appends its placeholder. A
   false return restores n_seq_map to seq_to_cnt.n, one per sequence. *)
AddOk == NewNSeq(Recs[cur]) < MaxNSeq
Add ==
    /\ group' = Append(group, cur)
    /\ hapSeq' = Refine(Recs[cur].v)
    /\ nseq' = NewNSeq(Recs[cur])
    /\ lowFreq' = Append(lowFreq, 0)

(* flush_compressed, then pc next: seq_coder3_flush hands out the coder's
   records, and initializes the coder unless it held none; low_freq moves
   to ready with the k-th placeholder replaced by the k-th coded record. A
   placeholder beyond the coded records exits. *)
Flushed(lf, g) ==
    [p \in 1..Len(lf) |->
        IF lf[p] # 0 THEN Item(lf[p], {})
        ELSE Item(g[Cardinality({q \in 1..p : lf[q] = 0})], Range(g))]
FlushThen(next) ==
    IF Placeholders > Len(group)
    THEN /\ out' = Append(out, Throw) /\ pc' = "throw"
         /\ UNCHANGED <<lowFreq, group, hapSeq, nseq, ready>>
    ELSE /\ ready' = ready \o Flushed(lowFreq, group)
         /\ lowFreq' = <<>>
         /\ group' = <<>>
         /\ hapSeq' = IF group = <<>> THEN hapSeq ELSE Fresh
         /\ nseq' = IF group = <<>> THEN nseq ELSE 1
         /\ pc' = next
         /\ UNCHANGED out

TypeOK ==
    /\ Recs \in [1..N -> Records]
    /\ i \in 1..(N + 1)
    /\ nblk \in Nat
    /\ blk \in Seq(1..N) /\ Len(blk) <= B + 1
    /\ bj \in 0..Len(blk)
    /\ cur \in 0..N
    /\ lastChrom \in 0..C
    /\ lowFreq \in Seq(0..N)
    /\ group \in Seq(1..N)
    /\ hapSeq \in [Haps -> Haps]
    /\ nseq \in Nat
    /\ dropped \in BOOLEAN
    /\ pc \in {"fill", "read", "add", "code", "retry", "ret", "idle", "throw"}

Init ==
    /\ Recs \in [1..N -> Records]
    /\ i = 1
    /\ nblk = 0
    /\ blk = <<>>
    /\ bj = 0
    /\ cur = 0
    /\ lastChrom = 0
    /\ lowFreq = <<>>
    /\ group = <<>>
    /\ hapSeq = Fresh
    /\ nseq = 1
    /\ ready = <<>>
    /\ pending = None
    /\ dropped = FALSE
    /\ pc = "fill"
    /\ out = <<>>
    /\ expected = Expected

Done == pc = "throw" \/ Len(out) = Len(expected)

(* The records fill's current block can take: RefIt.combine puts the
   first data line in front of BlockLineReader's first block. *)
Cap == IF nblk = 1 THEN B + 1 ELSE B

(* fill's loop with ready empty starts a block. *)
NewBlock ==
    /\ pc = "fill" /\ ready = <<>>
    /\ pc' = "read" /\ blk' = <<>> /\ bj' = 0 /\ nblk' = nblk + 1
    /\ UNCHANGED <<Recs, i, cur, lastChrom, lowFreq, group, hapSeq, nseq, ready, pending,
                   dropped, out, expected>>

(* The block's loop: block_reader_next gives a record. *)
Read ==
    /\ pc = "read" /\ i <= N /\ Len(blk) < Cap
    /\ blk' = Append(blk, i) /\ i' = i + 1
    /\ UNCHANGED <<Recs, nblk, bj, cur, lastChrom, lowFreq, group, hapSeq, nseq, ready,
                   pending, dropped, pc, out, expected>>

(* The block is full or block_reader_next gave NULL; an empty block
   flushes and ends fill. *)
ReadEnd ==
    /\ pc = "read" /\ (i > N \/ Len(blk) = Cap)
    /\ IF blk = <<>>
         THEN FlushThen("ret")
         ELSE /\ pc' = "add"
              /\ UNCHANGED <<lowFreq, group, hapSeq, nseq, ready, out>>
    /\ UNCHANGED <<Recs, i, nblk, blk, bj, cur, lastChrom, pending, dropped, expected>>

(* add_rec: the filter may drop the record; one on a new chromosome
   flushes the coder. *)
AddRec ==
    /\ pc = "add" /\ bj < Len(blk)
    /\ bj' = bj + 1
    /\ LET r == blk[bj + 1]
       IN  IF Recs[r].kind = "excl"
             THEN UNCHANGED <<cur, lastChrom, lowFreq, group, hapSeq, nseq, ready, pc, out>>
           ELSE /\ cur' = r /\ lastChrom' = Recs[r].chrom
                /\ IF lastChrom # 0 /\ Recs[r].chrom # lastChrom
                     THEN FlushThen("code")
                     ELSE /\ pc' = "code"
                          /\ UNCHANGED <<lowFreq, group, hapSeq, nseq, ready, out>>
    /\ UNCHANGED <<Recs, i, nblk, blk, pending, dropped, expected>>

(* The block is coded; fill's loop tests ready again. *)
AddEnd ==
    /\ pc = "add" /\ bj = Len(blk)
    /\ pc' = "fill"
    /\ UNCHANGED <<Recs, i, nblk, blk, bj, cur, lastChrom, lowFreq, group, hapSeq, nseq,
                   ready, pending, dropped, out, expected>>

(* apply_seq_coding false: the record goes to low_freq. True: it goes to
   the coder, and a full coder is flushed before the retry. *)
Code ==
    /\ pc = "code"
    /\ IF Recs[cur].kind = "low"
         THEN /\ lowFreq' = Append(lowFreq, cur) /\ pc' = "add"
              /\ UNCHANGED <<group, hapSeq, nseq, ready, out>>
       ELSE IF AddOk
         THEN /\ Add /\ pc' = "add"
              /\ UNCHANGED <<ready, out>>
       ELSE FlushThen("retry")
    /\ UNCHANGED <<Recs, i, nblk, blk, bj, cur, lastChrom, pending, dropped, expected>>

(* The add after the flush. A second false drops the record and keeps its
   placeholder. *)
Retry ==
    /\ pc = "retry"
    /\ IF AddOk
         THEN /\ Add /\ UNCHANGED dropped
         ELSE /\ lowFreq' = Append(lowFreq, 0) /\ dropped' = TRUE
              /\ UNCHANGED <<group, hapSeq, nseq>>
    /\ pc' = "add"
    /\ UNCHANGED <<Recs, i, nblk, blk, bj, cur, lastChrom, ready, pending, out, expected>>

(* fill returns, and ref_it_next returns the record it took. *)
Return ==
    /\ pc = "ret" \/ (pc = "fill" /\ ready # <<>>)
    /\ pc' = "idle"
    /\ out' = IF pending = None THEN out ELSE Append(out, pending)
    /\ pending' = None
    /\ UNCHANGED <<Recs, i, nblk, blk, bj, cur, lastChrom, lowFreq, group, hapSeq, nseq,
                   ready, dropped, expected>>

(* ref_it_next: NULL when ready is empty, else its first record, calling
   fill when that empties ready. *)
Take ==
    /\ pc = "idle" /\ ~Done
    /\ IF ready = <<>>
         THEN /\ out' = Append(out, End)
              /\ UNCHANGED <<ready, pending, pc>>
       ELSE IF Len(ready) = 1
         THEN /\ pending' = Head(ready) /\ ready' = <<>> /\ pc' = "fill"
              /\ UNCHANGED out
       ELSE /\ out' = Append(out, Head(ready)) /\ ready' = Tail(ready)
            /\ UNCHANGED <<pending, pc>>
    /\ UNCHANGED <<Recs, i, nblk, blk, bj, cur, lastChrom, lowFreq, group, hapSeq, nseq,
                   dropped, expected>>

Stop == Done /\ pc \in {"idle", "throw"} /\ UNCHANGED vars

Next == NewBlock \/ Read \/ ReadEnd \/ AddRec \/ AddEnd \/ Code \/ Retry \/ Return \/ Take \/ Stop

Spec == Init /\ [][Next]_vars
        /\ WF_vars(NewBlock \/ Read \/ ReadEnd \/ AddRec \/ AddEnd \/ Code \/ Retry \/ Return \/ Take)

-----------------------------------------------------------------------------
(* The caller sees what the reference gives: the kept records in file
   order, each coded record with its group, and then NULL, or as many of
   them as Java 5.4 returns before the flush that throws. *)
OutPrefix ==
    /\ Len(out) <= Len(expected)
    /\ \A k \in 1..Len(out) : out[k] = expected[k]

(* Each placeholder in low_freq has one coded record in the coder, but for
   the one a refused retry left. *)
PlaceholdersFilled == Placeholders = Len(group) + IF dropped THEN 1 ELSE 0

(* n_seq_map counts the distinct allele strings of the coder's records,
   one per sequence a haplotype carries, and stays below max_nseq. *)
CoderSeqs ==
    /\ nseq = NSeq(group)
    /\ nseq = Cardinality(Seqs)
    /\ nseq < MaxNSeq

(* An empty coder is in its initial state, which seq_coder3_flush relies
   on when it returns early. *)
EmptyCoderFresh == group = <<>> => hapSeq = Fresh /\ nseq = 1

(* A group never spans two chromosomes. *)
GroupOneChrom == \A j, k \in Range(group) : Recs[j].chrom = Recs[k].chrom

(* The caller sees every record and then NULL, or the exit. *)
Completes == <>(out = expected)
=============================================================================
