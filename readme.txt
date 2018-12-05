This archive contains two parts: 
(1) Some real infrared images with small targets and an implementation
(2) The small target detecton algorithm presented in our TIP paper [1].

If you use this code in your publications, please cite:
[1] Chenqiang Gao, Deyu Meng, Yi Yang, et al., "Infrared Patch-Image Model for Small Target Detection in a Single Image, " Image Processing, IEEE Transactions on, vol. 22, no. 12, pp. 4996-5009, 2013. 
@article{Gao2013IPI,
   author = {Gao, Chenqiang and Meng, Deyu and Yang, Yi and Wang, Yongtao and Zhou, Xiaofang and Hauptmann, Alex},
   title = {Infrared Patch-Image Model for Small Target Detection in a Single Image},
   journal = {Image Processing, IEEE Transactions on},
   volume = {22},
   number = {12},
   pages = {4996-5009},
   year = {2013}
}

If you use the test images in your publications, besides above reference, please cite the following reference, too:
[2] Gao Chenqiang, Zhang Tianqi,Li Qiang, ¡°Small infrared target detection using sparse ring representation,¡± IEEE Aerospace and Electronic Systems Magazine, vol. 27, no. 3, pp. 21-30, 2012. 
@article{Gao2012,
   author = {Chenqiang, Gao and Tianqi, Zhang and Qiang, Li},
   title = {Small infrared target detection using sparse ring representation},
   journal = {IEEE Aerospace and Electronic Systems Magazine},
   volume = {27},
   number = {3},
   pages = {21-30},
   year = {2012}
}

This code is just the core implementation. For real applicaton, we highly suggest to contain effective preprocessing and postprocessing steps for better performance.
Besides, a good parameter setting for the code discussed in our paper is also suggested to be carefully considered.

Please note that this code does not contain the segmentation step. If you want to get the locations of small targets, you need to further use some segmentation algorithm.

How to use this code?
--You just run the the file of main.m in matlab. Our matlab version for development is Matlab2013a. 

If you have any questions, please contact:
Author: Chenqiang Gao
Email: gaochenqiang@gmail.com or gaocq@cqupt.edu.cn
Copyright: Chongqing University of Posts and Telecommunications

* License: Our code is only available for non-commercial research use.