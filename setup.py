from setuptools import setup, find_packages
import os
from pathlib import Path

__version__ = ""
for line in open('asopipe/__init__.py'):
    if line.startswith('__version__'):
        exec(line.strip())

# setup pkgs to install
PACKAGES = find_packages()

# list for pre-requisite modules
# requirments.txt is preferred to pin specific version of modules
# ────────────────────────────────
# 2. requirements.txt 읽기
# ────────────────────────────────
def parse_requirements(fname="./config/requirements.txt"):
    req_path = Path(fname)
    if not req_path.is_file():        # 파일이 없으면 빈 리스트
        return []
    with req_path.open(encoding="utf-8") as fp:
        reqs = []
        for line in fp:
            line = line.strip()
            if line and not line.lstrip().startswith("#"):
                reqs.append(line)
        return reqs

REQUIRES = parse_requirements()

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files('asopipe')

setup(
    name='asopipe',
    version=__version__,
    description='package for common modules',
    author='dwkim',
    packages = PACKAGES,
    install_requires = REQUIRES,
    python_requires = '>=3.8',
    package_data={'':extra_files}
)