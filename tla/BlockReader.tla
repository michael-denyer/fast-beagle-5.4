----------------------------- MODULE BlockReader -----------------------------
(***************************************************************************)
(* The three-stage pipeline in src/vcf/block_reader.c. The reader thread   *)
(* pops a slot off the free stack, fills it with the next batch of lines   *)
(* (the batch after the last line has n = 0 and is the end sentinel) and   *)
(* appends it to the read FIFO; the parser thread pops read, parses, and   *)
(* appends to the full FIFO; the consumer thread, in block_reader_next,    *)
(* returns the batch it holds to the free stack and takes the head of      *)
(* full. The reader and parser return after the sentinel or when          *)
(* block_reader_close sets stop; the consumer's wait has no stop escape.   *)
(* Each step below is one critical section under the mutex, or the        *)
(* unlocked read_lines or parse_batch. All three threads wait on the one   *)
(* condition variable changed: a thread with pc "wait" is in its wait set, *)
(* and only a broadcast (Wake) or a spurious wakeup takes it out, so a     *)
(* lost wakeup shows up as a liveness violation.                           *)
(*                                                                         *)
(* Batches are numbered 1 .. NBatches; Sentinel = NBatches + 1 is the      *)
(* batch with n = 0. MayClose = FALSE checks the run with no close.        *)
(***************************************************************************)
EXTENDS Integers, Sequences, FiniteSets

CONSTANTS Slots, NBatches, MayClose

ASSUME Slots \in Nat \ {0} /\ NBatches \in Nat /\ MayClose \in BOOLEAN

SlotIds == 1 .. Slots
None == 0
Sentinel == NBatches + 1
Batches == 1 .. Sentinel

VARIABLES
    freeS,       \* free_slots, a stack: the top is the last element
    readQ,       \* read, a FIFO of slots holding unparsed batches
    fullQ,       \* full, a FIFO of slots holding parsed batches
    content,     \* per slot: the batch last read into it, or None
    nextBatch,   \* the batch the reader fills next
    stop,
    rpc,         \* the reader: "check", "wait", "fill" or "done"
    rslot,       \* the slot the reader is filling, or None
    ppc,         \* the parser: "check", "wait", "parse" or "done"
    pslot,       \* the slot the parser is parsing, or None
    cpc,         \* the consumer: "idle" (between calls), "wait", "check", "join" or "finished"
    cur,         \* the slot the consumer holds, or None
    taken        \* batches the consumer has taken, in order

vars == <<freeS, readQ, fullQ, content, nextBatch, stop, rpc, rslot, ppc, pslot, cpc, cur, taken>>

Range(s) == {s[i] : i \in DOMAIN s}
Held == {rslot, pslot, cur} \ {None}
InQueue == Range(freeS) \cup Range(readQ) \cup Range(fullQ)
SeenEnd == taken = Sentinel

(* pthread_cond_broadcast(&changed): every waiter re-checks its condition. *)
Wake(pc) == IF pc = "wait" THEN "check" ELSE pc

TypeOK ==
    /\ freeS \in Seq(SlotIds) /\ readQ \in Seq(SlotIds) /\ fullQ \in Seq(SlotIds)
    /\ content \in [SlotIds -> {None} \cup Batches]
    /\ nextBatch \in Batches
    /\ stop \in BOOLEAN
    /\ rpc \in {"check", "wait", "fill", "done"}
    /\ rslot \in {None} \cup SlotIds
    /\ ppc \in {"check", "wait", "parse", "done"}
    /\ pslot \in {None} \cup SlotIds
    /\ cpc \in {"idle", "wait", "check", "join", "finished"}
    /\ cur \in {None} \cup SlotIds
    /\ taken \in 0 .. Sentinel

Init ==
    /\ freeS = [i \in 1 .. Slots |-> i]
    /\ readQ = <<>>
    /\ fullQ = <<>>
    /\ content = [s \in SlotIds |-> None]
    /\ nextBatch = 1
    /\ stop = FALSE
    /\ rpc = "check" /\ rslot = None
    /\ ppc = "check" /\ pslot = None
    /\ cpc = "idle" /\ cur = None
    /\ taken = 0

(* read_batches: wait for a free slot unless stopped, then pop one. *)
ReaderCheck ==
    /\ rpc = "check"
    /\ IF stop
         THEN rpc' = "done" /\ UNCHANGED <<freeS, rslot>>
       ELSE IF freeS = <<>>
         THEN rpc' = "wait" /\ UNCHANGED <<freeS, rslot>>
       ELSE /\ rslot' = freeS[Len(freeS)]
            /\ freeS' = SubSeq(freeS, 1, Len(freeS) - 1)
            /\ rpc' = "fill"
    /\ UNCHANGED <<readQ, fullQ, content, nextBatch, stop, ppc, pslot, cpc, cur, taken>>

(* read_lines, then read_batches: append to read, broadcast, return after the sentinel. *)
ReaderPublish ==
    /\ rpc = "fill"
    /\ content' = [content EXCEPT ![rslot] = nextBatch]
    /\ readQ' = Append(readQ, rslot)
    /\ nextBatch' = IF nextBatch = Sentinel THEN Sentinel ELSE nextBatch + 1
    /\ rpc' = IF nextBatch = Sentinel THEN "done" ELSE "check"
    /\ rslot' = None
    /\ ppc' = Wake(ppc)
    /\ cpc' = Wake(cpc)
    /\ UNCHANGED <<freeS, fullQ, stop, pslot, cur, taken>>

(* parse_batches: wait for a read batch unless stopped, then pop it. *)
ParserCheck ==
    /\ ppc = "check"
    /\ IF stop
         THEN ppc' = "done" /\ UNCHANGED <<readQ, pslot>>
       ELSE IF readQ = <<>>
         THEN ppc' = "wait" /\ UNCHANGED <<readQ, pslot>>
       ELSE /\ pslot' = Head(readQ)
            /\ readQ' = Tail(readQ)
            /\ ppc' = "parse"
    /\ UNCHANGED <<freeS, fullQ, content, nextBatch, stop, rpc, rslot, cpc, cur, taken>>

(* parse_batch, then parse_batches: append to full, broadcast, return after the sentinel. *)
ParserPublish ==
    /\ ppc = "parse"
    /\ fullQ' = Append(fullQ, pslot)
    /\ ppc' = IF content[pslot] = Sentinel THEN "done" ELSE "check"
    /\ pslot' = None
    /\ rpc' = Wake(rpc)
    /\ cpc' = Wake(cpc)
    /\ UNCHANGED <<freeS, readQ, content, nextBatch, stop, rslot, cur, taken>>

(* block_reader_next's wait loop: wait while full is empty, else take its
   head. Nothing waits for full to shrink, so taking needs no broadcast. *)
TakeOrWait ==
    IF fullQ = <<>>
      THEN /\ cpc' = "wait"
           /\ cur' = None
           /\ UNCHANGED <<fullQ, taken>>
      ELSE /\ cur' = Head(fullQ)
           /\ fullQ' = Tail(fullQ)
           /\ taken' = taken + 1
           /\ cpc' = "idle"

(* block_reader_next with cur exhausted: return cur to free and broadcast,
   then take or wait. After the sentinel the call returns without locking,
   so it is not a step here. *)
ConsumerNext ==
    /\ cpc = "idle" /\ ~SeenEnd
    /\ IF cur = None
         THEN UNCHANGED <<freeS, rpc, ppc>>
         ELSE /\ freeS' = Append(freeS, cur)
              /\ rpc' = Wake(rpc)
              /\ ppc' = Wake(ppc)
    /\ TakeOrWait
    /\ UNCHANGED <<readQ, content, nextBatch, stop, rslot, pslot>>

(* The consumer's re-check after a wakeup. *)
ConsumerCheck ==
    /\ cpc = "check"
    /\ TakeOrWait
    /\ UNCHANGED <<freeS, readQ, content, nextBatch, stop, rpc, rslot, ppc, pslot>>

(* block_reader_close: set stop, broadcast, then join both threads. *)
Close ==
    /\ MayClose
    /\ cpc = "idle"
    /\ stop' = TRUE
    /\ cpc' = "join"
    /\ rpc' = Wake(rpc)
    /\ ppc' = Wake(ppc)
    /\ UNCHANGED <<freeS, readQ, fullQ, content, nextBatch, rslot, pslot, cur, taken>>

Join ==
    /\ cpc = "join"
    /\ rpc = "done" /\ ppc = "done"
    /\ cpc' = "finished"
    /\ UNCHANGED <<freeS, readQ, fullQ, content, nextBatch, stop, rpc, rslot, ppc, pslot, cur, taken>>

(* pthread_cond_wait may return with no broadcast. *)
SpuriousReader ==
    /\ rpc = "wait" /\ rpc' = "check"
    /\ UNCHANGED <<freeS, readQ, fullQ, content, nextBatch, stop, rslot, ppc, pslot, cpc, cur, taken>>

SpuriousParser ==
    /\ ppc = "wait" /\ ppc' = "check"
    /\ UNCHANGED <<freeS, readQ, fullQ, content, nextBatch, stop, rpc, rslot, pslot, cpc, cur, taken>>

SpuriousConsumer ==
    /\ cpc = "wait" /\ cpc' = "check"
    /\ UNCHANGED <<freeS, readQ, fullQ, content, nextBatch, stop, rpc, rslot, ppc, pslot, cur, taken>>

(* All three parties are finished: both threads returned, and the consumer
   either joined them or saw the sentinel and never closes. *)
AllDone ==
    /\ rpc = "done" /\ ppc = "done"
    /\ cpc = "finished" \/ (cpc = "idle" /\ SeenEnd)

Finished == AllDone /\ UNCHANGED vars

Next ==
    \/ ReaderCheck \/ ReaderPublish \/ SpuriousReader
    \/ ParserCheck \/ ParserPublish \/ SpuriousParser
    \/ ConsumerNext \/ ConsumerCheck \/ SpuriousConsumer
    \/ Close \/ Join
    \/ Finished

(* Threads run; the consumer keeps calling next until it closes or sees the
   sentinel. Close and spurious wakeups are not fair. *)
Fairness ==
    /\ WF_vars(ReaderCheck) /\ WF_vars(ReaderPublish)
    /\ WF_vars(ParserCheck) /\ WF_vars(ParserPublish)
    /\ WF_vars(ConsumerNext) /\ WF_vars(ConsumerCheck) /\ WF_vars(Join)

Spec == Init /\ [][Next]_vars /\ Fairness

-----------------------------------------------------------------------------
(* Every slot is in exactly one place: a queue, or held by one thread. *)
SlotsPartition ==
    /\ Len(freeS) + Len(readQ) + Len(fullQ) + Cardinality(Held) = Slots
    /\ InQueue \cup Held = SlotIds

(* The reader never fills a slot that is queued or held by another thread. *)
ReaderSlotExclusive ==
    rpc = "fill" => rslot \notin InQueue \cup {cur, pslot}

(* Batches reach the consumer in file order: full, then the parser's batch,
   then read hold the batches after the last one taken, in order. *)
Pipeline == fullQ \o (IF pslot = None THEN <<>> ELSE <<pslot>>) \o readQ
InOrder ==
    /\ \A i \in 1 .. Len(Pipeline) : content[Pipeline[i]] = taken + i
    /\ cur # None => content[cur] = taken

(* Once close is called, both threads return and the join completes. *)
CloseTerminates == [](cpc = "join" => <>(cpc = "finished"))

(* The consumer's wait in block_reader_next always ends. *)
WaitEnds == [](cpc = "wait" => <>(cpc = "idle"))

(* Every run ends with the sentinel seen or the reader closed. *)
Progress == <>(SeenEnd \/ cpc = "finished")

(* With MayClose = FALSE: the consumer sees the sentinel. *)
SeesSentinel == <>SeenEnd
=============================================================================
