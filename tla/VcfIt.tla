-------------------------------- MODULE VcfIt --------------------------------
(***************************************************************************)
(* The emission buffer of src/vcf/vcf_it.c: read_line (line 61),           *)
(* read_block (87), fill_buffer (123), vcf_it_next (139) and the           *)
(* fill_buffer in vcf_it_open (190). The iterator holds the next unread    *)
(* line (line, have_line). fill_buffer reads blocks of up to Batch lines   *)
(* until Batch records have passed the marker filter or the input ends,    *)
(* and the buffer holds the records that passed. read_block runs under     *)
(* util_try, so a failed read ends the block and is raised after the       *)
(* block's lines are parsed; a line that fails to parse is raised first.   *)
(* vcf_it_next takes the next buffered record and, when that empties the   *)
(* buffer, calls fill_buffer before it returns the record, so an error in  *)
(* that fill loses it. Lines 1..N are the data lines after the header;     *)
(* vcf_header_read exits when there are none, so N >= 1. Init chooses each *)
(* line's fate (Kind) and the line whose following read fails (FailAt, 0   *)
(* for none). Reading line 1 happens in vcf_header_read, outside util_try, *)
(* and is not modeled. Each step is one line read, one line parsed and     *)
(* filtered, or one call.                                                  *)
(***************************************************************************)
EXTENDS Naturals, Sequences, FiniteSets

CONSTANTS N, Batch

ASSUME N >= 1 /\ Batch >= 1

VARIABLES
    Kind,      \* per line: "ok", "excl" (fails the marker filter) or "bad"
    FailAt,    \* reading the line after this one fails; 0 for none
    pos,       \* the file line in it->line
    have,      \* it->have_line
    block,     \* it->lines[0..n_lines), as file line numbers
    j,         \* the lines of block that fill_buffer has parsed and filtered
    readErr,   \* the read_error of this block is not NULL
    buf,       \* it->buf[head..head + n), as file line numbers
    pending,   \* the record vcf_it_next returns once fill_buffer returns; 0 for none
    pc,        \* "fill" (the loop test), "read" (in read_block), "scan", "idle" or "exit"
    out        \* what vcf_it_next returned or raised, in order

vars == <<Kind, FailAt, pos, have, block, j, readErr, buf, pending, pc, out>>

Rec(i) == [k |-> "rec", line |-> i]
ParseErr(i) == [k |-> "parse", line |-> i]
ReadErr(i) == [k |-> "read", line |-> i]
End == [k |-> "end", line |-> 0]
Raised(o) == o.k \in {"parse", "read"}

-----------------------------------------------------------------------------
(* The reference: Java 5.4 VcfIt.fillEmissionBuffer and next, with its
   string buffer at the full 1024 lines (stringBufferSize is smaller only
   for a first data line longer than a sixteenth of the heap over 2048
   characters). The lines fall into fixed blocks of Batch, since every
   block but the last is full. A fill reads and parses whole blocks until Batch records have
   passed the filter or the lines run out, and the first error stops it.
   next() returns a record only after the fill it starts has returned.
   The error comes from the first block that holds a line that fails to
   parse or to be read past. fillStringBuffer reads the whole block before
   any of it is parsed, so a failed read in that block comes before its bad
   lines. *)
Blk(i) == (i - 1) \div Batch + 1
NBlk == Blk(N)
Lines(b) == {i \in 1..N : Blk(i) = b}
Min(S) == CHOOSE i \in S : \A k \in S : i <= k

Bad == {i \in 1..N : Kind[i] = "bad"}
(* The first line reading or parsing fails at; N + 1 for none. *)
ErrLine == Min(Bad \cup {i \in 1..N : i = FailAt} \cup {N + 1})
TheErr == IF FailAt # 0 /\ Blk(FailAt) = Blk(ErrLine) THEN ReadErr(FailAt) ELSE ParseErr(ErrLine)

Ok(a, s) == {i \in 1..N : a <= Blk(i) /\ Blk(i) <= s /\ Kind[i] = "ok"}

Seq1(S) == [m \in 1..Cardinality(S) |-> CHOOSE i \in S : Cardinality({k \in S : k <= i}) = m]
Recs(S) == [m \in 1..Cardinality(S) |-> Rec(Seq1(S)[m])]
Carry(c) == IF c = 0 THEN <<>> ELSE <<Rec(c)>>

(* What the caller sees from the fill that starts at block a, when the
   previous fill left record c (0 for none) to be returned after it. *)
RECURSIVE Ref(_, _)
Ref(a, c) ==
    IF a > NBlk THEN Carry(c) \o <<End>>
    ELSE LET s == Min({t \in a..NBlk : Cardinality(Ok(a, t)) >= Batch} \cup {NBlk})
             r == Recs(Ok(a, s))
         IN  IF ErrLine <= N /\ Blk(ErrLine) <= s THEN <<TheErr>>
             ELSE IF r = <<>> THEN Carry(c) \o <<End>>
             ELSE Carry(c) \o SubSeq(r, 1, Len(r) - 1) \o Ref(s + 1, r[Len(r)].line)

Expected == Ref(1, 0)

-----------------------------------------------------------------------------
TypeOK ==
    /\ Kind \in [1..N -> {"ok", "excl", "bad"}]
    /\ FailAt \in 0..N
    /\ pos \in 1..N + 1
    /\ have \in BOOLEAN
    /\ block \in Seq(1..N) /\ Len(block) <= Batch
    /\ j \in 0..Len(block)
    /\ readErr \in BOOLEAN
    /\ buf \in Seq(1..N)
    /\ pending \in 0..N
    /\ pc \in {"fill", "read", "scan", "idle", "exit"}

Init ==
    /\ Kind \in [1..N -> {"ok", "excl", "bad"}]
    /\ FailAt \in 0..N
    /\ pos = 1
    /\ have = TRUE
    /\ block = <<>>
    /\ j = 0
    /\ readErr = FALSE
    /\ buf = <<>>
    /\ pending = 0
    /\ pc = "fill"
    /\ out = <<>>

Done == pc = "exit" \/ (out # <<>> /\ out[Len(out)] = End)

(* fill_buffer's loop test: read another block, or return, and then
   vcf_it_next returns the record it took. *)
FillTest ==
    /\ pc = "fill"
    /\ IF have /\ Len(buf) < Batch
         THEN /\ pc' = "read" /\ block' = <<>> /\ j' = 0 /\ readErr' = FALSE
              /\ UNCHANGED <<out, pending>>
         ELSE /\ pc' = "idle"
              /\ out' = out \o Carry(pending)
              /\ pending' = 0
              /\ UNCHANGED <<block, j, readErr>>
    /\ UNCHANGED <<Kind, FailAt, pos, have, buf>>

(* One pass of the read_block loop: line moves into the block and
   read_line reads the one after it, which fails after line FailAt. The
   failure longjmps out of read_line to the util_try in fill_buffer. *)
ReadLine ==
    /\ pc = "read" /\ Len(block) < Batch /\ have
    /\ block' = Append(block, pos)
    /\ IF pos = FailAt
         THEN /\ readErr' = TRUE /\ pc' = "scan" /\ UNCHANGED <<pos, have>>
         ELSE /\ pos' = pos + 1 /\ have' = (pos + 1 <= N) /\ UNCHANGED <<readErr, pc>>
    /\ UNCHANGED <<Kind, FailAt, j, buf, pending, out>>

(* read_block returns: the block is full or the file ended. *)
ReadEnd ==
    /\ pc = "read" /\ (Len(block) = Batch \/ ~have)
    /\ pc' = "scan"
    /\ UNCHANGED <<Kind, FailAt, pos, have, block, j, readErr, buf, pending, out>>

(* fill_buffer raises the block's read error before it parses the block. *)
ReadFail ==
    /\ pc = "scan" /\ readErr
    /\ out' = Append(out, ReadErr(block[Len(block)])) /\ pc' = "exit"
    /\ UNCHANGED <<Kind, FailAt, pos, have, block, j, readErr, buf, pending>>

(* fill_buffer's loop over the parsed block: a parse error exits, and a
   record that passes the filter joins the buffer. *)
Scan ==
    /\ pc = "scan" /\ ~readErr /\ j < Len(block)
    /\ LET i == block[j + 1]
       IN  CASE Kind[i] = "bad" -> /\ out' = Append(out, ParseErr(i)) /\ pc' = "exit"
                                   /\ UNCHANGED buf
             [] Kind[i] = "ok"  -> /\ buf' = Append(buf, i) /\ UNCHANGED <<out, pc>>
             [] OTHER           -> UNCHANGED <<buf, out, pc>>
    /\ j' = j + 1
    /\ UNCHANGED <<Kind, FailAt, pos, have, block, readErr, pending>>

(* The block is used up, and the loop goes on. *)
ScanEnd ==
    /\ pc = "scan" /\ ~readErr /\ j = Len(block)
    /\ pc' = "fill"
    /\ UNCHANGED <<Kind, FailAt, pos, have, block, j, readErr, buf, pending, out>>

(* vcf_it_next: NULL when the buffer is empty, else its first record,
   calling fill_buffer when that empties the buffer. *)
Call ==
    /\ pc = "idle" /\ ~Done
    /\ IF buf = <<>>
         THEN /\ out' = Append(out, End) /\ UNCHANGED <<buf, pending, pc>>
       ELSE IF Len(buf) = 1
         THEN /\ pending' = Head(buf) /\ buf' = <<>> /\ pc' = "fill"
              /\ UNCHANGED out
       ELSE /\ out' = Append(out, Rec(Head(buf))) /\ buf' = Tail(buf)
            /\ UNCHANGED <<pending, pc>>
    /\ UNCHANGED <<Kind, FailAt, pos, have, block, j, readErr>>

Stop == Done /\ pc \in {"idle", "exit"} /\ UNCHANGED vars

Next == FillTest \/ ReadLine \/ ReadEnd \/ ReadFail \/ Scan \/ ScanEnd \/ Call \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(FillTest \/ ReadLine \/ ReadEnd \/ ReadFail \/ Scan \/ ScanEnd \/ Call)

-----------------------------------------------------------------------------
(* The caller sees what the reference gives for the same lines: the
   records of each completed fill, the last one only once the next fill
   returns, and then the end or the first error. *)
OutPrefix ==
    /\ Len(out) <= Len(Expected)
    /\ \A k \in 1..Len(out) : out[k] = Expected[k]

(* it->buf has room for 2 * BUFFER_SIZE records (vcf_it.c:189): a block
   starts only below BUFFER_SIZE and adds at most BUFFER_SIZE. *)
BufFits == Len(buf) <= 2 * Batch

(* The caller sees every record and then the end or the first error. *)
Completes == <>(out = Expected)
=============================================================================
