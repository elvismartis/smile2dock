import setuptools
import some_build_toolkit  # comes from the `some-build-toolkit` library

def get_version():
    version = some_build_toolkit.compute_version()
    return version

setuptools.setup(
    name="my-project",
    version=get_version(),
)
