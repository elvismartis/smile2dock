from setuptools import setup, find_packages
import pathlib
import io

HERE = pathlib.Path(__file__).parent

def read_file(name: str) -> str:
	p = HERE / name
	if p.exists():
		return p.read_text(encoding="utf-8")
	return ""

def parse_requirements() -> list:
	reqs = []
	text = read_file("requirements.txt")
	if not text:
		return reqs
	for line in text.splitlines():
		line = line.strip()
		if not line or line.startswith('#'):
			continue
		reqs.append(line)
	return reqs

long_description = read_file("README.md")
install_requires = parse_requirements()

setup(
	name="smile2dock",
	version="3.1.0",
	description="SMILES to 3D converter with optional GNN implicit-solvent minimization",
	long_description=long_description,
	long_description_content_type="text/markdown",
	author="elvismartis",
	author_email="elvis.afmartis@gmail.com",
	url="https://github.com/elvismartis/smile2dock",
	packages=find_packages(exclude=("tests", "examples", "obsolete_files", "update_patch")),
	py_modules=[
		"smile2dock_v3",
		"smile2dock_v3_GNNimplicitsolvent",
	],
	include_package_data=True,
	install_requires=install_requires,
	extras_require={
		"gnn": [
			"gnnimplicitsolvent @ git+https://github.com/fjclark/GNNImplicitSolvent.git",
			"torch>=2.0.0",
			"torch-geometric>=2.3.0",
			"openmm>=8.0.0",
			"openff-toolkit>=0.15.0",
		],
		"dev": ["pytest>=7.0.0", "flake8>=5.0.0"],
	},
	entry_points={
		"console_scripts": [
			"smile2dock=smile2dock_v3:main",
			"smile2dock-gnn=smile2dock_v3_GNNimplicitsolvent:main",
		]
	},
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	python_requires=">=3.8",
)

