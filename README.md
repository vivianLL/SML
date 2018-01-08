# SML
class model for Supervised Multiclass Labeling(SML) which is a text annotation algorithm for images.
corel5k中部分图像基于SML算法的类模型

getDCT.m 实现单个图像的YBR空间转换，分块抽取DCT系数，并对系数进行降维  

gauss.m 为main.m的子函数，计算高斯分布密度函数  

cmnpdf.m 为gmm.m的子函数，计算多元正态概率密度函数  
covfixer2.m 为gmm.m的子函数，使协方差矩阵强制转换为有效的协方差矩阵  
getargs.m 为gmm.m的子函数，对函数的默认参数进行处理  

gmm.m 实现将每一幅图像聚类成8个混合高斯分布  

main.m 为主文件，结合上述文件中函数功能，读取文件夹中所有的图像，实现将多个高斯分布聚类成64个混合高斯分布  

最终聚类的混合高斯分布参数存储在model.mat文件中。  

mat文件中包含64个高斯分布的μ值，sigma值,和π值。  


说明：一个图像块的数据是192维的，为了使得提高程序的运行速率本程序中将数据压缩到30维。 

具体原理详见博客：
有监督的多类标注（SML）的原理及matlab实现：http://blog.csdn.net/vivian_ll/article/details/79004473
使用EM算法估计GMM参数的原理及matlab实现：http://blog.csdn.net/vivian_ll/article/details/78793293
