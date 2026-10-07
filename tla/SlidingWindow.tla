---------------------------- MODULE SlidingWindow ----------------------------
(***************************************************************************)
(* The read-ahead hand-over in src/vcf/sliding_window.c: hand_over,        *)
(* read_windows, read_ahead, sliding_window_next and sliding_window_close. *)
(* The reader thread produces windows 1..NWindows and then NULL (0 here),  *)
(* and hands each over through ahead and ready, waiting until the caller   *)
(* clears ready or sets stop. Producing window k with k in FailAt may      *)
(* instead fail (util_try catches util_exit): the reader sets error and    *)
(* ready and returns. The caller takes windows with sliding_window_next    *)
(* and may call sliding_window_close at any point, which sets stop and     *)
(* joins the reader. Each step below is one critical section under the     *)
(* mutex, or the unlocked production of a window, as in BlockReader.tla.   *)
(* Both threads wait on the one condition variable changed: a thread with  *)
(* pc "wait" is in its wait set, and only a broadcast (Wake) or a spurious *)
(* wakeup takes it out, so a lost wakeup shows up as a liveness violation. *)
(***************************************************************************)
EXTENDS Integers, Sequences

CONSTANTS NWindows, FailAt

ASSUME NWindows \in Nat \ {0} /\ FailAt \subseteq 1 .. NWindows

Null == 0
Windows == 1 .. NWindows

VARIABLES
    ready,     \* sw->ready
    stop,      \* sw->stop
    error,     \* sw->error != NULL
    ahead,     \* sw->ahead: Null or a window
    rpc,       \* the reader: "produce", "hand", "fail", "wait", "check" or "done"
    w,         \* the reader's local w in read_windows: Null or a window
    produced,  \* windows made so far (sw->window_index)
    failedAt,  \* the window whose production failed, or Null
    cpc,       \* the consumer: "idle" (between calls), "wait", "check", "join", "closed" or "exited"
    got        \* what sliding_window_next returned, in order: windows or Null

vars == <<ready, stop, error, ahead, rpc, w, produced, failedAt, cpc, got>>

Range(s) == {s[i] : i \in DOMAIN s}

(* pthread_cond_broadcast(&changed): every waiter re-checks its condition. *)
Wake(pc) == IF pc = "wait" THEN "check" ELSE pc

TypeOK ==
    /\ ready \in BOOLEAN /\ stop \in BOOLEAN /\ error \in BOOLEAN
    /\ ahead \in {Null} \cup Windows
    /\ rpc \in {"produce", "hand", "fail", "wait", "check", "done"}
    /\ w \in {Null} \cup Windows
    /\ produced \in 0 .. NWindows
    /\ failedAt \in {Null} \cup Windows
    /\ cpc \in {"idle", "wait", "check", "join", "closed", "exited"}
    /\ got \in Seq({Null} \cup Windows) /\ Len(got) <= NWindows + 1

Init ==
    /\ ready = FALSE /\ stop = FALSE /\ error = FALSE
    /\ ahead = Null
    /\ rpc = "produce" /\ w = Null /\ produced = 0 /\ failedAt = Null
    /\ cpc = "idle" /\ got = <<>>

-----------------------------------------------------------------------------
(* Reader thread *)

(* read_windows: next_*_window makes window produced + 1, or fails inside
   it, or w stays NULL once the last window has been read. *)
Produce ==
    /\ rpc = "produce"
    /\ IF produced = NWindows
         THEN w' = Null /\ rpc' = "hand" /\ UNCHANGED <<produced, failedAt>>
         ELSE \/ w' = produced + 1 /\ produced' = produced + 1 /\ rpc' = "hand" /\ UNCHANGED failedAt
              \/ produced + 1 \in FailAt /\ failedAt' = produced + 1 /\ rpc' = "fail" /\ UNCHANGED <<w, produced>>
    /\ UNCHANGED <<ready, stop, error, ahead, cpc, got>>

(* hand_over's wait loop, then read_windows' return test. *)
AwaitTaken(isReady) ==
    rpc' = IF isReady /\ ~stop THEN "wait"
           ELSE IF stop \/ w = Null THEN "done" ELSE "produce"

(* hand_over up to its first wait: publish w and broadcast. *)
HandOver ==
    /\ rpc = "hand"
    /\ ahead' = w /\ ready' = TRUE
    /\ cpc' = Wake(cpc)
    /\ AwaitTaken(TRUE)
    /\ UNCHANGED <<stop, error, w, produced, failedAt, got>>

(* hand_over's re-check after a wakeup. *)
ReaderCheck ==
    /\ rpc = "check"
    /\ AwaitTaken(ready)
    /\ UNCHANGED <<ready, stop, error, ahead, w, produced, failedAt, cpc, got>>

(* read_ahead after util_try returns an error. *)
Fail ==
    /\ rpc = "fail"
    /\ error' = TRUE /\ ready' = TRUE
    /\ cpc' = Wake(cpc)
    /\ rpc' = "done"
    /\ UNCHANGED <<stop, ahead, w, produced, failedAt, got>>

-----------------------------------------------------------------------------
(* Consumer: the main thread *)

(* sliding_window_next under the mutex: wait, or exit on error, or take ahead
   (clearing ready and waking the reader), or return NULL leaving ready set.
   The caller may keep calling next after NULL and gets NULL each time; got
   records the first NULL only, which keeps the state space finite. *)
TakeOrWait ==
    IF ~ready
      THEN cpc' = "wait" /\ UNCHANGED <<ready, ahead, rpc, got>>
    ELSE IF error
      THEN cpc' = "exited" /\ UNCHANGED <<ready, ahead, rpc, got>>
    ELSE IF ahead # Null
      THEN /\ got' = Append(got, ahead)
           /\ ahead' = Null /\ ready' = FALSE
           /\ rpc' = Wake(rpc)
           /\ cpc' = "idle"
    ELSE /\ got' = IF Null \in Range(got) THEN got ELSE Append(got, Null)
         /\ cpc' = "idle"
         /\ UNCHANGED <<ready, ahead, rpc>>

ConsumerNext ==
    /\ cpc \in {"idle", "check"}
    /\ TakeOrWait
    /\ UNCHANGED <<stop, error, w, produced, failedAt>>

(* sliding_window_close: set stop, broadcast, then join the reader. *)
Close ==
    /\ cpc = "idle"
    /\ stop' = TRUE
    /\ rpc' = Wake(rpc)
    /\ cpc' = "join"
    /\ UNCHANGED <<ready, error, ahead, w, produced, failedAt, got>>

Join ==
    /\ cpc = "join" /\ rpc = "done"
    /\ cpc' = "closed"
    /\ UNCHANGED <<ready, stop, error, ahead, rpc, w, produced, failedAt, got>>

(* pthread_cond_wait may return with no broadcast. *)
SpuriousReader ==
    /\ rpc = "wait" /\ rpc' = "check"
    /\ UNCHANGED <<ready, stop, error, ahead, w, produced, failedAt, cpc, got>>

SpuriousConsumer ==
    /\ cpc = "wait" /\ cpc' = "check"
    /\ UNCHANGED <<ready, stop, error, ahead, rpc, w, produced, failedAt, got>>

Finished == cpc \in {"closed", "exited"} /\ UNCHANGED vars

Next ==
    \/ Produce \/ HandOver \/ ReaderCheck \/ Fail \/ SpuriousReader
    \/ ConsumerNext \/ Close \/ Join \/ SpuriousConsumer
    \/ Finished

(* Threads run; a started next runs to its return. Whether the caller calls
   next or close, and spurious wakeups, are not fair. *)
Fairness ==
    /\ WF_vars(Produce) /\ WF_vars(HandOver) /\ WF_vars(ReaderCheck) /\ WF_vars(Fail)
    /\ WF_vars(cpc = "check" /\ ConsumerNext) /\ WF_vars(Join)

Spec == Init /\ [][Next]_vars /\ Fairness

-----------------------------------------------------------------------------
(* The consumer receives windows 1, 2, 3, ... in order, none skipped or
   repeated, and NULL only after window NWindows. *)
InOrder ==
    \A i \in DOMAIN got : got[i] = i \/ (got[i] = Null /\ i > NWindows)

(* Once a window is in the consumer's hands, ahead never points at it, so
   sliding_window_close's window_free(sw->ahead) is never a double free. *)
AheadNotGiven == ahead # Null => ahead \notin Range(got)

(* After join, every window made is either with the consumer or in ahead,
   where close frees it. *)
NoLeak == cpc = "closed" => \A k \in 1 .. produced : k \in Range(got) \/ k = ahead

(* The reader fails producing window k+1 only after the consumer has taken
   windows 1..k, and the consumer takes nothing more once error is set. *)
ErrorAfterPrefix == error => Len(got) = failedAt - 1 /\ Null \notin Range(got)

(* The reader never leaves a window behind when it fails. *)
ErrorNoAhead == error => ahead = Null

(* Liveness *)
CloseTerminates == [](cpc = "join" => <>(cpc = "closed"))
NextReturns == [](cpc \in {"wait", "check"} => <>(cpc \in {"idle", "exited"}))
ExitOnlyOnError == [](cpc = "exited" => error /\ Len(got) = failedAt - 1)
=============================================================================
