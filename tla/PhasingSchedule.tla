--------------------------- MODULE PhasingSchedule ---------------------------
(***************************************************************************)
(* The iteration schedule of one window: phase_window and                  *)
(* phase_stage1_variants in src/main/main.c, phase_data_init,              *)
(* phase_data_increment_it and phase_data_advance_to_first_phasing_it in   *)
(* src/phase/phase_data.c, the EM gating and PBWT direction of             *)
(* phase_ls_run_stage1 in src/phase/phase_ls.c, and the burn-in flag that  *)
(* phase_baum1_init (src/phase/phase_baum1.c) sets from pd->it, which      *)
(* decides whether update_phase shrinks the unphased lists. par.c allows   *)
(* burnin >= 1 and iterations >= 1, as Java's Par does. Init chooses       *)
(* whether the input is already phased, whether every marker is a stage-1  *)
(* marker, and par->em; each run's swap rate is chosen when it ends: "low" *)
(* (at most 1%), "high", or "nan" (0 / 0 when no sample had an unphased    *)
(* het, for which <= is false), or "unread" once the run reached burnin,   *)
(* where && never reaches the comparison. Each step is one run of          *)
(* phase_ls_run_stage1 with the increment after it, the swap-rate test     *)
(* that follows, or one path choice in phase_window.                       *)
(***************************************************************************)
EXTENDS Integers, Sequences

CONSTANTS Burnin, Iterations

ASSUME Burnin \in Nat \ {0} /\ Iterations \in Nat \ {0}

VARIABLES
    phased,     \* targ_is_phased(w)
    allStage1,  \* fpd->n_stage1 == fpd->n_markers
    emOn,       \* p->em
    it,         \* pd->it
    pc,         \* "start", "loop", "check", "stage2", "write" or "done"
    runs,       \* each phase_ls_run_stage1: its it, pb->burnin, EM step, direction
    rates,      \* each run's swap rate
    stage2It,   \* pd->it when phase_ls_run_stage2 starts; -1 if it does not
    path        \* "none", "phased", "stage1" or "stage2"

vars == <<phased, allStage1, emOn, it, pc, runs, rates, stage2It, path>>

NIts == Burnin + Iterations
Rates == {"low", "high", "nan"}

(* swap <= max_burnin_swap_rate *)
AtMostOnePercent(r) == r = "low"

(* while (pd->it < n_its) *)
LoopCond == it < NIts

Its == [i \in 1 .. Len(runs) |-> runs[i].it]

(* The runs a window must make when burn-in ends after k runs: burn-in
   iterations 0 .. k-1, then every phasing iteration once, in order. *)
Schedule(k) == [i \in 1 .. k |-> i - 1] \o [i \in 1 .. Iterations |-> Burnin + i - 1]

IsPrefix(s, t) == Len(s) <= Len(t) /\ \A i \in 1 .. Len(s) : s[i] = t[i]

TypeOK ==
    /\ phased \in BOOLEAN /\ allStage1 \in BOOLEAN /\ emOn \in BOOLEAN
    /\ it \in 0 .. NIts
    /\ pc \in {"start", "loop", "check", "stage2", "write", "done"}
    /\ runs \in Seq([it : 0 .. NIts - 1, burnin : BOOLEAN, em : {"init", "update", "none"}, bwd : BOOLEAN])
    /\ rates \in Seq(Rates \cup {"unread"})
    /\ stage2It \in -1 .. NIts
    /\ path \in {"none", "phased", "stage1", "stage2"}

Init ==
    /\ phased \in BOOLEAN /\ allStage1 \in BOOLEAN /\ emOn \in BOOLEAN
    /\ it = 0
    /\ pc = "start"
    /\ runs = <<>> /\ rates = <<>>
    /\ stage2It = -1 /\ path = "none"

(* phase_window: the input is already phased, or stage 1 runs. *)
Begin ==
    /\ pc = "start"
    /\ IF phased THEN pc' = "write" /\ path' = "phased" ELSE pc' = "loop" /\ UNCHANGED path
    /\ UNCHANGED <<phased, allStage1, emOn, it, runs, rates, stage2It>>

(* One pass of the while loop: phase_ls_run_stage1 at it, with EM
   initialised at it 0 and updated during burn-in, the backward PBWT on even
   iterations, and each phase_baum1's burnin flag; then
   phase_data_increment_it. *)
Run ==
    /\ pc = "loop" /\ LoopCond
    /\ runs' = Append(runs, [it |-> it, burnin |-> it < Burnin,
                             em |-> IF ~emOn THEN "none" ELSE IF it = 0 THEN "init"
                                    ELSE IF it < Burnin THEN "update" ELSE "none",
                             bwd |-> it % 2 = 0])
    /\ it' = it + 1
    /\ pc' = "check"
    /\ UNCHANGED <<phased, allStage1, emOn, rates, stage2It, path>>

(* The swap rate of the run and the early end of burn-in. *)
Check ==
    /\ pc = "check"
    /\ \E r \in IF it < Burnin THEN Rates ELSE {"unread"} :
         /\ rates' = Append(rates, r)
         /\ IF it < Burnin /\ AtMostOnePercent(r)
              THEN it' = Burnin
              ELSE UNCHANGED it
    /\ pc' = "loop"
    /\ UNCHANGED <<phased, allStage1, emOn, runs, stage2It, path>>

(* The loop ends; phase_window picks the stage-1 or the stage-2 haplotypes. *)
LoopEnd ==
    /\ pc = "loop" /\ ~LoopCond
    /\ IF allStage1 THEN pc' = "write" /\ path' = "stage1" ELSE pc' = "stage2" /\ path' = "stage2"
    /\ UNCHANGED <<phased, allStage1, emOn, it, runs, rates, stage2It>>

(* phase_ls_run_stage2, whose pbwt_ibs_data_init reads pd->it. *)
Stage2 ==
    /\ pc = "stage2"
    /\ stage2It' = it /\ pc' = "write"
    /\ UNCHANGED <<phased, allStage1, emOn, it, runs, rates, path>>

(* The window is written. *)
Write ==
    /\ pc = "write" /\ pc' = "done"
    /\ UNCHANGED <<phased, allStage1, emOn, it, runs, rates, stage2It, path>>

Stop == pc = "done" /\ UNCHANGED vars

Next == Begin \/ Run \/ Check \/ LoopEnd \/ Stage2 \/ Write \/ Stop

Spec == Init /\ [][Next]_vars /\ WF_vars(Begin \/ Run \/ Check \/ LoopEnd \/ Stage2 \/ Write)

-----------------------------------------------------------------------------
(* At most burnin + iterations runs. *)
RunBound == Len(runs) <= NIts

(* The runs so far are a prefix of a schedule, and a finished stage 1 ran a
   whole one: the burn-in runs come first, and every phasing iteration runs
   once, in order. *)
FollowsSchedule ==
    /\ \E k \in 1 .. Burnin : IsPrefix(Its, Schedule(k))
    /\ (pc \in {"stage2", "write", "done"} /\ ~phased) => \E k \in 1 .. Burnin : Its = Schedule(k)

(* Burn-in leaves early only after a rate of at most 1%, and always does. *)
EarlyExitOnLowRate ==
    \A i \in 1 .. Len(runs) - 1 :
        (runs[i + 1].it # runs[i].it + 1 <=> (rates[i] = "low" /\ runs[i].it + 1 < Burnin))

(* update_phase leaves the unphased lists alone in every burn-in run and
   shrinks them in every phasing run. *)
BurninRuns == \A i \in 1 .. Len(runs) : runs[i].burnin <=> runs[i].it < Burnin

(* A finished stage 1 ends with iterations runs that shrink the unphased
   lists, the last run among them. *)
LastRunsPhase ==
    (pc = "done" /\ ~phased) =>
        /\ ~runs[Len(runs)].burnin
        /\ \A i \in 1 .. Len(runs) : ~runs[i].burnin <=> i > Len(runs) - Iterations

(* With em=true, EM initialises the parameters in the first run and updates
   them only in the other burn-in runs; with em=false it never runs. *)
EmInBurnin ==
    \A i \in 1 .. Len(runs) :
        /\ (runs[i].em = "init") <=> (emOn /\ i = 1)
        /\ (runs[i].em = "update") <=> (emOn /\ i > 1 /\ runs[i].it < Burnin)

(* Stage 2 starts after the last increment, where pbwt_ibs_data_init takes
   the stage-2 candidate count, and only for a window with markers outside
   stage 1. *)
Stage2AfterStage1 ==
    /\ stage2It \in {-1, NIts}
    /\ path = "stage2" => ~phased /\ ~allStage1
    /\ pc = "done" => (stage2It = NIts <=> ~phased /\ ~allStage1)

(* A phased input runs nothing. *)
PhasedRunsNothing == phased => runs = <<>>

Terminates == <>(pc = "done")
=============================================================================
