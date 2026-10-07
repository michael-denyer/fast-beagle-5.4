----------------------------- MODULE BgenCleanup -----------------------------
(***************************************************************************)
(* The partial-output cleanup in src/bgen/bgen_files.c. bgen_files_open    *)
(* creates a BGEN member and records it, bgen_files_complete marks the     *)
(* members complete, and the atexit handler remove_partial removes the     *)
(* recorded members of a run that did not complete. Each is one critical   *)
(* section under the file's lock, so each is one step here. exit() can     *)
(* start on any thread at any point (tla/FatalExit.tla) while the others   *)
(* keep running, so an open or a complete can come after it.               *)
(***************************************************************************)
CONSTANTS Members

VARIABLES
    proc,       \* "running", "exiting" (inside exit()), "exit1" or "exit0"
    created,    \* members bgen_files_open created
    completed,  \* bgen_files_complete ran
    cleaned,    \* remove_partial ran
    removed     \* members remove_partial removed

vars == <<proc, created, completed, cleaned, removed>>

TypeOK ==
    /\ proc \in {"running", "exiting", "exit1", "exit0"}
    /\ created \subseteq Members
    /\ completed \in BOOLEAN
    /\ cleaned \in BOOLEAN
    /\ removed \subseteq Members

Init ==
    /\ proc = "running"
    /\ created = {}
    /\ completed = FALSE
    /\ cleaned = FALSE
    /\ removed = {}

(* The finished flag of bgen_files.c. *)
Finished == completed \/ cleaned
Live == proc \in {"running", "exiting"}

(* bgen_files_open: creates nothing once finished. *)
Open(m) ==
    /\ Live /\ m \notin created /\ ~Finished
    /\ created' = created \cup {m}
    /\ UNCHANGED <<proc, completed, cleaned, removed>>

(* bgen_files_complete, from window_writer_close. Once finished is set the
   call changes nothing. *)
Complete ==
    /\ Live /\ ~Finished
    /\ completed' = TRUE
    /\ UNCHANGED <<proc, created, cleaned, removed>>

(* A thread calls exit(). *)
Exit ==
    /\ proc = "running"
    /\ proc' = "exiting"
    /\ UNCHANGED <<created, completed, cleaned, removed>>

(* exit() runs remove_partial while the other threads still run. *)
Cleanup ==
    /\ proc = "exiting" /\ ~cleaned
    /\ removed' = IF completed THEN {} ELSE created
    /\ cleaned' = TRUE
    /\ UNCHANGED <<proc, created, completed>>

(* exit() then ends the process with status 1. *)
End ==
    /\ proc = "exiting" /\ cleaned
    /\ proc' = "exit1"
    /\ UNCHANGED <<created, completed, cleaned, removed>>

(* main returns 0 after window_writer_close. *)
Return ==
    /\ proc = "running" /\ completed
    /\ proc' = "exit0"
    /\ UNCHANGED <<created, completed, cleaned, removed>>

Ended == ~Live /\ UNCHANGED vars

Next == (\E m \in Members : Open(m)) \/ Complete \/ Exit \/ Cleanup \/ End \/ Return \/ Ended

Spec == Init /\ [][Next]_vars /\ WF_vars(Cleanup) /\ WF_vars(End)

-----------------------------------------------------------------------------
(* Once remove_partial ran, every member on disk is complete: it removed the
   partial ones and nothing is created after it. *)
NoPartialAfterCleanup == cleaned => created \subseteq removed \/ completed

(* bgen_files_complete before exit keeps the members. *)
CompletedSurvive == completed => removed = {}

(* exit() always finishes. *)
Terminates == [](proc = "exiting" => <>(proc = "exit1"))
=============================================================================
