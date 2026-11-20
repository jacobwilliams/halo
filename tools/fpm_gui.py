"""Run fpm for halo with the C++ settings that its GUI dependencies need.

Run inside the pixi environment, from the halo root, e.g.:

    pixi run gui-build                               (= python tools/fpm_gui.py build)
    pixi run gui                                     (= python tools/fpm_gui.py run halo_gui)
    python tools/fpm_gui.py run halo_gui -- examples/example_sparse.json   (the -- is optional)
    python tools/fpm_gui.py run halo_solver -- examples/example_sparse.json

vkframes and qt-fortran each have a C++ layer that fpm compiles (VulkanSceneGraph
and Qt). Their own helper scripts supply the compiler and flags through fpm's
environment variables (FPM_CXX, FPM_CXXFLAGS, FPM_LDFLAGS):

  * vkframes' tools/fpm_vsg.py: VSG's flags, the compiled shaders, and
    VKFRAMES_DATA_PATH (its star catalog and textures);
  * qt-fortran's tools/fpm_qt.py: Qt's flags (it keeps the ones already set),
    and runs fpm.

Since they are fpm dependencies, both scripts are in build/dependencies (fetched
by `fpm update` if missing).

The fpm command is run with the release profile and OpenMP, as in run.sh
(without qr_mumps). Extra fpm flags can be given in HALO_FPM_FLAGS.
"""

import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # halo
DEPS = os.path.join(ROOT, "build", "dependencies")
USAGE = __doc__


def main(argv):
    if not argv or argv[0] not in ("build", "run"):
        sys.exit(USAGE)
    os.chdir(ROOT)
    if not (os.path.isdir(os.path.join(DEPS, "vkframes")) and os.path.isdir(os.path.join(DEPS, "qt-fortran"))):
        subprocess.check_call(["fpm", "update"])

    sys.path.insert(0, os.path.join(DEPS, "vkframes", "tools"))
    import fpm_vsg  # noqa: E402  (vkframes' build settings)

    env_settings, _ = fpm_vsg.settings()
    env = dict(os.environ, **env_settings)

    flags = ["--profile", "release", "--flag", ("-fopenmp " + os.environ.get("HALO_FPM_FLAGS", "")).strip()]
    if argv[0] == "build":
        fpm = ["build", *flags, *argv[1:]]
    else:
        # run <target> [--] [args...]: everything after the target is the program's
        # (the `--` is optional: `pixi run gui -- args` drops it)
        rest = argv[1:]
        target = []
        if rest and rest[0] != "--":
            target = [rest[0]]
            rest = rest[1:]
        if rest and rest[0] == "--":
            rest = rest[1:]
        fpm = ["run", *target, *flags, "--", *rest]

    qt_tool = os.path.join(DEPS, "qt-fortran", "tools", "fpm_qt.py")
    return subprocess.call([sys.executable, qt_tool, *fpm], env=env)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
