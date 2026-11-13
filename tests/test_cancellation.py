"""Signal and worker cancellation qualification; no fixed millisecond guarantee."""

import json
import os
import signal
import subprocess
import sys
import sysconfig
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest
from cave_closure import CancellationToken, IncompleteResultError, prepare_json
from generate_graphs import generate_graph


def test_free_threaded_import_does_not_enable_gil():
    if sysconfig.get_config_var("Py_GIL_DISABLED"):
        assert not sys._is_gil_enabled()


def test_worker_cancellation_during_owned_computation():
    problem = generate_graph(
        seed=713,
        stations=3000,
        extra_edges=0,
        topology="chain",
        controls_per_component=2,
    )
    prepared = prepare_json(json.dumps(problem))
    token = CancellationToken()
    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(
            prepared.solve, reduce=False, cancellation=token, max_iterations=100000
        )
        time.sleep(0.03)
        token.cancel()
        try:
            result = future.result(timeout=30)
        except IncompleteResultError:
            return
    assert result.status == "Cancelled"
    with pytest.raises(RuntimeError):
        result.require_converged()


@pytest.mark.skipif(os.name == "nt", reason="POSIX SIGINT qualification")
def test_main_thread_keyboard_interrupt_is_preserved():
    script = """
import json, sys
sys.path.insert(0, sys.argv[1])
from generate_graphs import generate_graph
from cave_closure import prepare_json
prepared = prepare_json(json.dumps(generate_graph(seed=713, stations=3000, extra_edges=0, topology="chain")))
print("ready", flush=True)
try:
    prepared.solve(reduce=False, max_iterations=100000)
except KeyboardInterrupt:
    print("interrupted", flush=True)
    raise SystemExit(0)
raise SystemExit(2)
"""
    process = subprocess.Popen(
        [sys.executable, "-c", script, str(Path(__file__).parent)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        assert process.stdout.readline().strip() == "ready"
        time.sleep(0.1)
        process.send_signal(signal.SIGINT)
        output, error = process.communicate(timeout=30)
        assert process.returncode == 0, error
        assert "interrupted" in output
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
