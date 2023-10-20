from  setuptools import setup, find_packages

with open("PYPI.rst", "r") as f:
    long_description = f.read()

setup(
    name='stereosite',
    version='1.0.5',
    author='LiuXing',
    author_email='liuxing2@genomics.cn',
    description=('Analysis spatial transcriptomics data'),
    long_description=long_description,
    license='MIT License',
    keywords='spatial cell interaction intensity',
    url="https://github.com/STOmics/StereoSiTE",

    packages=['stereosite'], #需要打包的目录列表

    include_package_data=True,
    platforms='any',
    #需要安装的依赖包
    install_requires = [
        'anndata>=0.8.0',
        'scanpy>=1.9.1',
        'squidpy>=1.1.2',
        'decoupler>=1.4.0',
        'pydeseq2>=0.3.6',
        'networkx>=3.1'
    ],
)