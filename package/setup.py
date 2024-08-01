from setuptools import setup
import re
def extract_version(file_path):
    """
    Extracts the version string from a given file.

    Parameters:
    file_path (str): The path to the file from which to extract the version.

    Returns:
    str: The extracted version string, or None if the version string is not found.
    """
    version_pattern = r"^__version__ = ['\"]([^'\"]*)['\"]"

    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(version_pattern, line)
            if match:
                return match.group(1)

    raise RuntimeError("Version could not be obtained. Check the format of the __version__.py file.")


version = extract_version("src/sails/__version__.py")

setup(
    cmake_source_dir=".",
    version=version
)