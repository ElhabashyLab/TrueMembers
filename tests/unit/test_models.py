import pytest
from dataclasses import FrozenInstanceError

from truemembers.models import (
    DetectionMethod,
    Hit,
    IterationSummary,
    JobResult,
    JobStatus,
    JobSummary,
    RunResult,
    StopReason,
)


def make_hit(**overrides) -> Hit:
    values = {
        "uniprot_id": "P00533",
        "region_start": 1,
        "region_end": 100,
        "detection_method": DetectionMethod.JACKHMMER,
    }
    values.update(overrides)
    return Hit(**values)


class TestDetectionMethod:
    def test_members_and_values(self):
        assert DetectionMethod.JACKHMMER.value == "JACKHMMER"
        assert DetectionMethod.FOLDSEEK.value == "FOLDSEEK"
        assert DetectionMethod.JACKHMMER_AND_FOLDSEEK.value == "JACKHMMER_AND_FOLDSEEK"

    def test_merged_same_method_is_identity(self):
        assert DetectionMethod.JACKHMMER.merged(DetectionMethod.JACKHMMER) is DetectionMethod.JACKHMMER
        assert DetectionMethod.FOLDSEEK.merged(DetectionMethod.FOLDSEEK) is DetectionMethod.FOLDSEEK

    def test_merged_different_methods_combines(self):
        assert (
            DetectionMethod.JACKHMMER.merged(DetectionMethod.FOLDSEEK)
            is DetectionMethod.JACKHMMER_AND_FOLDSEEK
        )
        assert (
            DetectionMethod.FOLDSEEK.merged(DetectionMethod.JACKHMMER)
            is DetectionMethod.JACKHMMER_AND_FOLDSEEK
        )

    def test_merged_combined_absorbs_other_methods(self):
        for method in DetectionMethod:
            assert (
                DetectionMethod.JACKHMMER_AND_FOLDSEEK.merged(method)
                is DetectionMethod.JACKHMMER_AND_FOLDSEEK
            )


class TestHit:
    def test_defaults(self):
        hit = make_hit()
        assert hit.unique_id
        assert hit.bitscore is None
        assert hit.tm_score_to_target is None
        assert hit.ali_coverage is None
        assert hit.domains == {}
        assert hit.is_all_member is False
        assert hit.is_strict_member is False
        assert hit.first_seen_iteration == 0
        assert hit.structure_path is None
        assert hit.sequence_path is None

    def test_unique_ids_are_unique(self):
        assert make_hit().unique_id != make_hit().unique_id

    def test_domains_are_not_shared_between_instances(self):
        first = make_hit()
        second = make_hit()
        first.domains["PF00069"] = (1, 10)
        assert second.domains == {}

    def test_region_property(self):
        assert make_hit(region_start=507, region_end=569).region == (507, 569)

    def test_is_frozen(self):
        hit = make_hit()
        with pytest.raises(FrozenInstanceError):
            hit.uniprot_id = "P04637"

    @pytest.mark.parametrize("region_start,region_end", [(0, 10), (-5, 10), (10, 5), (10, 1)])
    def test_invalid_regions_raise(self, region_start, region_end):
        with pytest.raises(ValueError):
            make_hit(region_start=region_start, region_end=region_end)

    def test_empty_uniprot_id_raises(self):
        with pytest.raises(ValueError):
            make_hit(uniprot_id="")


class TestJobResult:
    def test_defaults(self):
        result = JobResult(seed=make_hit())
        assert result.hits == []
        assert result.sequence_hits == 0
        assert result.structure_hits == 0
        assert result.deduped_hits == 0
        assert result.downloads == 0
        assert result.structural_pass == 0
        assert result.failed is False

    def test_failed_when_error_recorded(self):
        assert JobResult(seed=make_hit(), error="boom").failed is True

    def test_hits_lists_are_not_shared(self):
        first = JobResult(seed=make_hit())
        second = JobResult(seed=make_hit())
        first.hits.append(make_hit())
        assert second.hits == []


class TestSummaries:
    def test_job_status_matches_seed_use_states(self):
        assert [status.value for status in JobStatus] == [
            "PENDING",
            "RUNNING",
            "COMPLETED",
            "FAILED",
            "DROPPED",
        ]

    def test_stop_reasons(self):
        assert [reason.value for reason in StopReason] == [
            "CONVERGED",
            "MAX_ITERATIONS",
            "SEEDS_EXHAUSTED",
            "FAILED",
        ]

    def test_job_summary_defaults(self):
        summary = JobSummary(
            job_id="job-1",
            iteration=1,
            seed_id="seed-1",
            seed_uniprot_id="Q5SK37",
            seed_region_start=507,
            seed_region_end=569,
        )
        assert summary.status is JobStatus.COMPLETED
        assert summary.attempts == 1
        assert summary.new_strict_members == 0
        assert summary.duration_s == 0.0
        assert summary.error is None

    def test_iteration_summary_duration(self):
        summary = IterationSummary(
            iteration=1,
            seeds_submitted=2,
            jobs_completed=1,
            jobs_failed=1,
            new_strict_members=3,
            total_strict_members=4,
            growth=0.5,
            started_at="2026-10-06T15:00:00",
            finished_at="2026-10-06T15:02:30",
        )
        assert summary.duration_s == 150.0

    @pytest.mark.parametrize(
        "started_at,finished_at",
        [(None, None), ("2026-10-06T15:00:00", None), (None, "2026-10-06T15:02:30")],
    )
    def test_iteration_summary_duration_none_without_timestamps(self, started_at, finished_at):
        summary = IterationSummary(
            iteration=1,
            seeds_submitted=1,
            jobs_completed=1,
            jobs_failed=0,
            new_strict_members=0,
            total_strict_members=0,
            started_at=started_at,
            finished_at=finished_at,
        )
        assert summary.duration_s is None

    def test_run_result_defaults_and_duration(self):
        result = RunResult(
            target_uid="Q5SK37",
            target_region_start=507,
            target_region_end=569,
            stop_reason=StopReason.CONVERGED,
            run_dir="/tmp/run",
            started_at="2026-10-06T15:00:00",
            finished_at="2026-10-06T15:05:00",
        )
        assert result.iterations == []
        assert result.jobs == []
        assert result.total_all_members == 0
        assert result.total_strict_members == 0
        assert result.config_snapshot == {}
        assert result.duration_s == 300.0

    def test_run_result_lists_are_not_shared(self):
        first = RunResult(
            target_uid="Q5SK37",
            target_region_start=507,
            target_region_end=569,
            stop_reason=StopReason.FAILED,
            run_dir="/tmp/first",
        )
        second = RunResult(
            target_uid="Q5SK37",
            target_region_start=507,
            target_region_end=569,
            stop_reason=StopReason.FAILED,
            run_dir="/tmp/second",
        )
        first.iterations.append(
            IterationSummary(
                iteration=1,
                seeds_submitted=0,
                jobs_completed=0,
                jobs_failed=0,
                new_strict_members=0,
                total_strict_members=0,
            )
        )
        assert second.iterations == []
