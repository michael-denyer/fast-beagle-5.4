--------------------------- MODULE ParallelOrdered ---------------------------
(***************************************************************************)
(* parallel_ordered in src/blbutil/parallel.c. Workers claim items in      *)
(* increasing order, but only while next < consumed + Window, and mark     *)
(* slot item % Window built. The calling thread consumes items in order,   *)
(* clearing each slot. Each step below is one critical section under the   *)
(* mutex, or the unlocked build or consume. A condition-variable wait adds *)
(* the thread to a wait set, and only a signal, a broadcast or a spurious  *)
(* wakeup takes it out, so a lost wakeup shows up as a liveness violation  *)
(* (a spurious wakeup is always possible, so never as a TLC deadlock).     *)
(***************************************************************************)
EXTENDS Integers, FiniteSets

CONSTANTS Workers, NItems, Window

ASSUME NItems \in Nat /\ Window \in Nat \ {0} /\ Workers # {}

Slots == 0 .. Window - 1
None == -1

VARIABLES
    next,        \* the next item to claim
    consumed,    \* items consumed, in item order
    slot,        \* per slot: the built, unconsumed item in it, or None
    wpc,         \* per worker: "check", "wait", "build" or "done"
    witem,       \* per worker: the item it is building
    roomWait,    \* workers waiting on room_cv
    cpc,         \* the consumer: "check", "wait", "consume", "join" or "finished"
    citem,       \* the consumer's loop item
    builtWait    \* the consumer is waiting on built_cv

vars == <<next, consumed, slot, wpc, witem, roomWait, cpc, citem, builtWait>>

TypeOK ==
    /\ next \in 0 .. NItems
    /\ consumed \in 0 .. NItems
    /\ slot \in [Slots -> {None} \cup (0 .. NItems - 1)]
    /\ wpc \in [Workers -> {"check", "wait", "build", "done"}]
    /\ witem \in [Workers -> {None} \cup (0 .. NItems - 1)]
    /\ roomWait \subseteq Workers
    /\ cpc \in {"check", "wait", "consume", "join", "finished"}
    /\ citem \in 0 .. NItems
    /\ builtWait \in BOOLEAN

Init ==
    /\ next = 0
    /\ consumed = 0
    /\ slot = [s \in Slots |-> None]
    /\ wpc = [w \in Workers |-> "check"]
    /\ witem = [w \in Workers |-> None]
    /\ roomWait = {}
    /\ cpc = "check"
    /\ citem = 0
    /\ builtWait = FALSE

(* run_ordered: the while test, then claim or finish. *)
WorkerCheck(w) ==
    /\ wpc[w] = "check"
    /\ IF next < NItems /\ next >= consumed + Window
         THEN /\ wpc' = [wpc EXCEPT ![w] = "wait"]
              /\ roomWait' = roomWait \cup {w}
              /\ UNCHANGED <<next, witem>>
         ELSE IF next >= NItems
         THEN /\ wpc' = [wpc EXCEPT ![w] = "done"]
              /\ UNCHANGED <<next, witem, roomWait>>
         ELSE /\ witem' = [witem EXCEPT ![w] = next]
              /\ next' = next + 1
              /\ wpc' = [wpc EXCEPT ![w] = "build"]
              /\ UNCHANGED roomWait
    /\ UNCHANGED <<consumed, slot, cpc, citem, builtWait>>

(* build() done: mark the slot, and signal built_cv if it is the consumer's item. *)
WorkerBuilt(w) ==
    /\ wpc[w] = "build"
    /\ slot' = [slot EXCEPT ![witem[w] % Window] = witem[w]]
    /\ IF witem[w] = consumed /\ builtWait
         THEN /\ builtWait' = FALSE
              /\ cpc' = "check"
         ELSE UNCHANGED <<builtWait, cpc>>
    /\ wpc' = [wpc EXCEPT ![w] = "check"]
    /\ witem' = [witem EXCEPT ![w] = None]
    /\ UNCHANGED <<next, consumed, roomWait, citem>>

(* The consumer's loop: wait for its item's slot, or leave to join. *)
ConsumerCheck ==
    /\ cpc = "check"
    /\ IF citem >= NItems
         THEN /\ cpc' = "join"
              /\ UNCHANGED builtWait
         ELSE IF slot[citem % Window] = None
         THEN /\ cpc' = "wait"
              /\ builtWait' = TRUE
         ELSE /\ cpc' = "consume"
              /\ UNCHANGED builtWait
    /\ UNCHANGED <<next, consumed, slot, wpc, witem, roomWait, citem>>

(* consume() done: clear the slot, count the item, broadcast room_cv. *)
Consume ==
    /\ cpc = "consume"
    /\ slot' = [slot EXCEPT ![citem % Window] = None]
    /\ consumed' = citem + 1
    /\ citem' = citem + 1
    /\ wpc' = [w \in Workers |-> IF w \in roomWait THEN "check" ELSE wpc[w]]
    /\ roomWait' = {}
    /\ cpc' = "check"
    /\ UNCHANGED <<next, witem, builtWait>>

Join ==
    /\ cpc = "join"
    /\ \A w \in Workers : wpc[w] = "done"
    /\ cpc' = "finished"
    /\ UNCHANGED <<next, consumed, slot, wpc, witem, roomWait, citem, builtWait>>

(* pthread_cond_wait may return with no signal. *)
SpuriousWorker(w) ==
    /\ w \in roomWait
    /\ roomWait' = roomWait \ {w}
    /\ wpc' = [wpc EXCEPT ![w] = "check"]
    /\ UNCHANGED <<next, consumed, slot, witem, cpc, citem, builtWait>>

SpuriousConsumer ==
    /\ builtWait
    /\ builtWait' = FALSE
    /\ cpc' = "check"
    /\ UNCHANGED <<next, consumed, slot, wpc, witem, roomWait, citem>>

Finished == cpc = "finished" /\ UNCHANGED vars

Next ==
    \/ \E w \in Workers : WorkerCheck(w) \/ WorkerBuilt(w) \/ SpuriousWorker(w)
    \/ ConsumerCheck \/ Consume \/ Join \/ SpuriousConsumer
    \/ Finished

Fairness ==
    /\ \A w \in Workers : WF_vars(WorkerCheck(w)) /\ WF_vars(WorkerBuilt(w))
    /\ WF_vars(ConsumerCheck) /\ WF_vars(Consume) /\ WF_vars(Join)

Spec == Init /\ [][Next]_vars /\ Fairness

-----------------------------------------------------------------------------
(* At most Window items are claimed but not consumed, so slots never clash. *)
WindowBound == next - consumed <= Window

(* A worker never overwrites a slot whose item has not been consumed. *)
NoOverwrite ==
    \A w \in Workers : wpc[w] = "build" => slot[witem[w] % Window] = None

(* The consumer consumes exactly its own item. *)
ConsumesOwnItem == cpc = "consume" => slot[citem % Window] = citem

(* The consumer takes items one at a time in order. *)
InOrder == consumed = citem

(* Every run ends with every item consumed and every worker joined. *)
Terminates == <>(cpc = "finished")
AllConsumed == [](cpc = "finished" => consumed = NItems)
=============================================================================
