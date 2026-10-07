------------------------------ MODULE FatalExit ------------------------------
(***************************************************************************)
(* The fatal-error lifecycle of util_exit and util_oom in                  *)
(* src/blbutil/utilities.c across threads. A thread inside util_try that   *)
(* calls util_exit longjmps back to its try frame with the message stored; *)
(* the consumer of that frame raises it later (sliding_window_next,        *)
(* block_reader_next, vcf_it_next). A thread outside util_try that calls   *)
(* util_exit, or any thread that calls util_oom, tests the exiting flag:   *)
(* set, it pauses forever; clear, it calls exit(), which ends the process. *)
(*                                                                         *)
(* A longjmp out of a critical section would leave the mutex locked for    *)
(* good, so the code raises input errors outside its locks. Under the      *)
(* chrom_ids lock only an allocation can fail, and util_oom does not       *)
(* unwind: the thread exits or pauses with the lock held, and the process  *)
(* ends without it.                                                        *)
(*                                                                         *)
(* Threads: main, the read-ahead reader (read_windows runs under util_try  *)
(* on it), and parallel_for workers that enter util_try per parse task.    *)
(***************************************************************************)
CONSTANTS Workers

ASSUME Workers # {}

Main == "main"
Reader == "reader"
Threads == {Main, Reader} \cup Workers
None == "none"

VARIABLES
    proc,       \* "running", "exiting" (inside exit()), "exit1" or "exit0"
    exiting,    \* the atomic_flag in util_exit
    pc,         \* per thread: "run", "wait" (on the mutex), "paused", "exit" or "done"
    inTry,      \* per thread: inside a util_try frame
    owner,      \* the thread holding the chrom_ids mutex, or None
    deferred    \* per thread: "none", "pending", "raised" or "discarded"

vars == <<proc, exiting, pc, inTry, owner, deferred>>

TypeOK ==
    /\ proc \in {"running", "exiting", "exit1", "exit0"}
    /\ exiting \in BOOLEAN
    /\ pc \in [Threads -> {"run", "wait", "paused", "exit", "done"}]
    /\ inTry \in [Threads -> BOOLEAN]
    /\ owner \in Threads \cup {None}
    /\ deferred \in [Threads -> {"none", "pending", "raised", "discarded"}]

Init ==
    /\ proc = "running"
    /\ exiting = FALSE
    /\ pc = [t \in Threads |-> "run"]
    /\ inTry = [t \in Threads |-> t = Reader]
    /\ owner = None
    /\ deferred = [t \in Threads |-> "none"]

Live == proc \in {"running", "exiting"}
(* The thread is running and not inside the chrom_ids critical section. *)
Unlocked(t) == Live /\ pc[t] = "run" /\ owner # t

(* chrom_ids_index and chrom_ids_id: lock, use the table, unlock.
   pthread_mutex_lock takes a free mutex, or blocks on a held one. *)
Lock(t) ==
    /\ Unlocked(t)
    /\ IF owner = None
         THEN owner' = t /\ UNCHANGED pc
         ELSE pc' = [pc EXCEPT ![t] = "wait"] /\ UNCHANGED owner
    /\ UNCHANGED <<proc, exiting, inTry, deferred>>

(* A blocked thread gets the mutex once it is free. *)
Acquire(t) ==
    /\ Live /\ pc[t] = "wait" /\ owner = None
    /\ owner' = t
    /\ pc' = [pc EXCEPT ![t] = "run"]
    /\ UNCHANGED <<proc, exiting, inTry, deferred>>

Unlock(t) ==
    /\ Live /\ pc[t] = "run" /\ owner = t
    /\ owner' = None
    /\ UNCHANGED <<proc, exiting, pc, inTry, deferred>>

(* parallel_for parse tasks run under util_try (parse_task in block_reader.c
   and vcf_it.c). The reader is under util_try for its whole life
   (read_ahead in sliding_window.c). *)
EnterTry(w) ==
    /\ Unlocked(w) /\ ~inTry[w]
    /\ inTry' = [inTry EXCEPT ![w] = TRUE]
    /\ UNCHANGED <<proc, exiting, pc, owner, deferred>>

LeaveTry(w) ==
    /\ Unlocked(w) /\ inTry[w]
    /\ inTry' = [inTry EXCEPT ![w] = FALSE]
    /\ UNCHANGED <<proc, exiting, pc, owner, deferred>>

(* util_exit for an input error inside util_try: store the message, longjmp
   to the try frame. The reader's frame is read_ahead, which parks the error
   and returns. A worker's frame is parse_task, which keeps the error for the
   consumer and goes on to the next item. *)
FailInTry(t) ==
    /\ Unlocked(t) /\ inTry[t]
    /\ deferred' = [deferred EXCEPT ![t] = "pending"]
    /\ inTry' = [inTry EXCEPT ![t] = FALSE]
    /\ pc' = [pc EXCEPT ![t] = IF t = Reader THEN "done" ELSE "run"]
    /\ UNCHANGED <<proc, exiting, owner>>

(* The tail of util_exit: the second exiting thread pauses forever; the
   first calls exit(). A lock the thread holds stays held. *)
ExitProcess(t) ==
    /\ exiting' = TRUE
    /\ IF exiting
         THEN pc' = [pc EXCEPT ![t] = "paused"] /\ UNCHANGED proc
         ELSE pc' = [pc EXCEPT ![t] = "exit"] /\ proc' = "exiting"

(* util_oom for a failed allocation, inside or outside util_try and with or
   without the lock, or util_exit for an input error outside util_try. *)
Fatal(t) ==
    /\ Live /\ pc[t] = "run"
    /\ ExitProcess(t)
    /\ UNCHANGED <<inTry, owner, deferred>>

(* exit() runs to the end and the process ends with status 1. *)
RunExit(t) ==
    /\ proc = "exiting" /\ pc[t] = "exit"
    /\ proc' = "exit1"
    /\ pc' = [pc EXCEPT ![t] = "done"]
    /\ UNCHANGED <<exiting, inTry, owner, deferred>>

(* The reader finishes its windows and returns normally. *)
ReaderDone ==
    /\ Unlocked(Reader)
    /\ pc' = [pc EXCEPT ![Reader] = "done"]
    /\ inTry' = [inTry EXCEPT ![Reader] = FALSE]
    /\ UNCHANGED <<proc, exiting, owner, deferred>>

(* A worker thread returns to parallel_for's join. *)
WorkerDone(w) ==
    /\ Unlocked(w) /\ ~inTry[w]
    /\ pc' = [pc EXCEPT ![w] = "done"]
    /\ UNCHANGED <<proc, exiting, inTry, owner, deferred>>

(* The reader takes a parsed block after parallel_for joined its workers and
   raises the first stored parse error (block_reader_next, vcf_it_next),
   still under its own util_try. *)
ConsumeWorker(w) ==
    /\ Unlocked(Reader) /\ inTry[Reader]
    /\ deferred[w] = "pending"
    /\ \A v \in Workers : pc[v] = "done"
    /\ deferred' = [deferred EXCEPT ![w] = "raised", ![Reader] = "pending"]
    /\ inTry' = [inTry EXCEPT ![Reader] = FALSE]
    /\ pc' = [pc EXCEPT ![Reader] = "done"]
    /\ UNCHANGED <<proc, exiting, owner>>

(* sliding_window_next raises the reader's parked error on the main thread,
   outside any util_try. *)
MainTake ==
    /\ Unlocked(Main)
    /\ deferred[Reader] = "pending"
    /\ deferred' = [deferred EXCEPT ![Reader] = "raised"]
    /\ ExitProcess(Main)
    /\ UNCHANGED <<inTry, owner>>

(* After the last window: sliding_window_close joins the reader and frees
   parse errors of blocks never taken (block_reader_close), main returns 0. *)
MainClose ==
    /\ Unlocked(Main)
    /\ pc[Reader] = "done" /\ deferred[Reader] # "pending"
    /\ \A w \in Workers : pc[w] = "done"
    /\ deferred' = [t \in Threads |-> IF deferred[t] = "pending" THEN "discarded" ELSE deferred[t]]
    /\ proc' = "exit0"
    /\ pc' = [pc EXCEPT ![Main] = "done"]
    /\ UNCHANGED <<exiting, inTry, owner>>

Ended == ~Live /\ UNCHANGED vars

Next ==
    \/ \E t \in Threads :
        Lock(t) \/ Acquire(t) \/ Unlock(t) \/ FailInTry(t) \/ Fatal(t) \/ RunExit(t)
    \/ \E w \in Workers : EnterTry(w) \/ LeaveTry(w) \/ WorkerDone(w) \/ ConsumeWorker(w)
    \/ ReaderDone \/ MainTake \/ MainClose \/ Ended

(* A critical section ends and exit() completes. The main thread gets the
   mutex however often the other threads take it, and calls
   sliding_window_next again however often it uses the lock in between. *)
Fairness ==
    /\ \A t \in Threads : WF_vars(Unlock(t)) /\ WF_vars(RunExit(t))
    /\ SF_vars(Acquire(Main)) /\ SF_vars(MainTake)

Spec == Init /\ [][Next]_vars /\ Fairness

-----------------------------------------------------------------------------
(* Once a thread reaches the tail of util_exit, the process ends. *)
Terminates == [](exiting => <>(proc = "exit1"))

(* The lock is released, or the process ends. *)
LockReleased == [](owner # None => <>(owner = None \/ ~Live))

(* The main thread never blocks on the lock for as long as the process lives. *)
MainNotBlocked == [](pc[Main] = "wait" => <>(pc[Main] # "wait" \/ ~Live))

(* The reader's deferred error is raised, unless another exit beats it. *)
DeferredRaised ==
    deferred[Reader] = "pending" ~> (deferred[Reader] = "raised" \/ proc = "exit1")
=============================================================================
