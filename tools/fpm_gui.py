"""Build and run the HALO GUI (the fpm package in gui/).

Run inside the pixi environment, from anywhere, e.g.:

    pixi run gui-build                         (= python tools/fpm_gui.py build)
    pixi run gui                               (= python tools/fpm_gui.py run)
    pixi run gui -- examples/example.json --run
    python tools/fpm_gui.py run [--] [halo_gui arguments...]

The solver itself (the halo root package) does not need any of this: it builds
and runs with plain fpm (see the README). The GUI is a separate package so that
the solver doesn't depend on Qt or Vulkan.

vkframes and qt-fortran each have a C++ layer that fpm compiles (VulkanSceneGraph
and Qt). Their own helper scripts supply the compiler and flags through fpm's
environment variables (FPM_CXX, FPM_CXXFLAGS, FPM_LDFLAGS):

  * vkframes' tools/fpm_vsg.py: VSG's flags, the compiled shaders, and
    VKFRAMES_DATA_PATH (its star catalog and textures);
  * qt-fortran's tools/fpm_qt.py: Qt's flags (it keeps the ones already set),
    and runs fpm.

Both are in gui/build/dependencies (fetched by `fpm update` if missing).

The GUI is built with the release profile and OpenMP, as in run.sh (without
qr_mumps); extra fpm flags can be given in HALO_FPM_FLAGS. It is run with
`--dir <halo root>`, since the config files' paths (data/, kernel/, ...) are
relative to it.
"""

import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # halo
GUI = os.path.join(ROOT, "gui")                                       # the GUI package
DEPS = os.path.join(GUI, "build", "dependencies")


def fpm_qt(env, command, extra):
    """Run fpm in the GUI package through qt-fortran's fpm_qt.py (which adds Qt's
    settings, and on Windows puts the DLLs' folders on PATH for `run`)."""
    flags = ["--profile", "release", "--flag", ("-fopenmp " + os.environ.get("HALO_FPM_FLAGS", "")).strip()]
    qt_tool = os.path.join(DEPS, "qt-fortran", "tools", "fpm_qt.py")
    return subprocess.call([sys.executable, qt_tool, command, *flags, *extra], env=env, cwd=GUI)


def main(argv):
    if not argv or argv[0] not in ("build", "run"):
        sys.exit(__doc__)
    if not (os.path.isdir(os.path.join(DEPS, "vkframes")) and os.path.isdir(os.path.join(DEPS, "qt-fortran"))):
        subprocess.check_call(["fpm", "update"], cwd=GUI)

    os.chdir(GUI)                                 # (fpm_vsg: the project's build dir is ./build)
    sys.path.insert(0, os.path.join(DEPS, "vkframes", "tools"))
    import fpm_vsg  # noqa: E402  (vkframes' build settings)

    env_settings, dll_dir = fpm_vsg.settings()
    env = dict(os.environ, **env_settings)
    if dll_dir:                                   # Windows: vkframes' DLL
        env["PATH"] = dll_dir + os.pathsep + env.get("PATH", "")

    if argv[0] == "build":
        return fpm_qt(env, "build", argv[1:])

    # run [--] [args...]  (pixi drops the --): run from the halo root (--dir)
    args = argv[1:]
    if args and args[0] == "--":
        args = args[1:]
    return fpm_qt(env, "run", ["halo_gui", "--", "--dir", ROOT, *args])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
