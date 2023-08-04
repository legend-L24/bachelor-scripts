诸位好，这部分内容是关于如何根据实验谱图解析MOF结构的代码（神经网络部分仍然存在问题），gasa_simulate.ipynb是如何拟合谱图。
首先，经过CNN文件夹中的代码构建神经网络，给出MOF结构可能的空间群
然后，在rough_build 文件夹中，利用贝叶斯优化晶胞参数，以及用decryst_new/src中基于c的代码，./solve + {指定结构}.txt优化金属原子坐标
接下来，在refine_code中精修金属原子位置
最后，在construct_mof中，自动识别金属原子位置并插入配体，生成接出的MOF的cif结构

