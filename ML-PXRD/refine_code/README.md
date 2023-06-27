Hello, the code use Bayesian Optimization to replace people to refine the cifs of MOFs structures. The package is based on https://www.nature.com/articles/s41524-020-0330-9 but some corrections are done since the difference between MOFs and mineral
代码运行案例
可在Wanggroup服务器/home/jyb/yifanhou/refine_code文件夹下面  
conda activate yfh  
python main.py  
重要的输入目录:  
可在main.py文件中直接调试  
存放cif文件的路径：/home/jyb/yifanhou/examples/1/cif  
存放目标谱图的目录：/home/jyb/yifanhou/examples/1/spec.csv  
存放仪器参数的目录：/home/jyb/yifanhou/examples/1/inst_xry.prm  

软件使用步骤：  
（1）激活环境，将多个可能的cif和一个谱图放入指定路径  
（2）运行python main.py，会输出 job_+{cif文件名}文件夹，里面储存着所有精修输出文件，Error_+{cif文件名}文本文件, 里面每一行是每一轮的最终误差。比如job_1.cif，Error_1.cif  
（3）运行python search.py，生成selected_and_refined文件夹，并挑选出最低误差的结果放入，如果search.py出现bug，请直接查找Error_+{cif文件名}文件，并根据序号查找job_+{cif文件名}文件夹中结果。  
可能的三种结果：  
损失函数高于某个阙值（search.py中设为20%），不显示  
损失函数低于阙值，如果能查找到坐标信息，输出结构谱图  
损失函数低于阙值，无法查找坐标信息，这是因为方法原理是选取其中一小段谱图修到最优并推广到整段谱图，该情况是属于推广失败，也算是精修失败。  
