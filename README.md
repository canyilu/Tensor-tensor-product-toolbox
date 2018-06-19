## Tensor-Tensor Product Toolbox

### Introduction

The tensor-tensor product (t-product) <a class="footnote-reference" href="#id2" id="id1">[1]</a> is a natural generalization of matrix multiplication. Based on t-product, many operations on matrix can be extended to tensor cases, including tensor SVD (see an illustration in the figure below), tensor spectral norm, tensor nuclear norm <a class="footnote-reference" href="#id2" id="id1">[2]</a> and many others. The linear algebraic structure of tensors are similar to the matrix cases. We develop a Matlab toolbox to implement several basic operations on tensors based on t-product.

<p align="center"> 
<img src="https://github.com/canyilu/tproduct/blob/master/tsvd.JPG">
</p>

### List of Functions

The table below gives the list of functions implemented in our toolbox. The detailed definitions of these tensor concepts, operations and tensor factorizations are given at <a href="../publications/2018-software-tproduct.pdf" class="textlink" target="_blank">https://canyilu.github.io/publications/2018-software-tproduct.pdf</a>. 

<p align="center"> 
<img src="https://github.com/canyilu/tproduct/blob/master/tab_funs.JPG">
</p>

Note that we only focus on 3 way tensor in this toolbox. We will develop the same functions for p-way tensor in the near future. We will also provide the python version soon.

Simply run the following routine to test all the above functions:
```matlab
test.m
```

### Citing

<p>In citing this toolbox in your papers, please use the following references:</p>

<div class="highlight-none"><div class="highlight"><pre>
C. Lu. Tensor-Tensor Product Toolbox. Carnegie Mellon University, June 2018. https://github.com/canyilu/tproduct.
C. Lu, J. Feng, Y. Chen, W. Liu, Z. Lin, and S. Yan. Tensor robust principal component analysis with a new tensor
nuclear norm. arXiv preprint arXiv:1804.03728, 2018.
</pre></div>

<p>The corresponding BiBTeX citations are given below:</p>
<div class="highlight-none"><div class="highlight"><pre>
@manual{lu2018tproduct,
  author       = {Lu, Canyi},
  title        = {Tensor-Tensor Product Toolbox},
  organization = {Carnegie Mellon University},
  month        = {June},
  year         = {2018},
  note         = {\url{https://github.com/canyilu/tproduct}}
}
@article{lu2018tensor,
  author       = {Lu, Canyi and Feng, Jiashi and Chen, Yudong and Liu, Wei and Lin, Zhouchen and Yan, Shuicheng},
  title        = {Tensor Robust Principal Component Analysis with A New Tensor Nuclear Norm},
  journal      = {arXiv preprint arXiv:1804.03728},
  year         = {2018}
}
</pre></div>
  
### Related Toolboxes
The t-product toolbox has been applied in our works about tensor roubst PCA <a class="footnote-reference" href="#id2" id="id1">[2,3]</a>, low-rank tensor completion and low-rank tensor recovery from Gaussian measurements <a class="footnote-reference" href="#id2" id="id1">[4]</a>. Some more models are included in LibADMM toolbox <a class="footnote-reference" href="#id2" id="id1">[5]</a>.
<ul>
  <li> <a href="https://github.com/canyilu/Tensor-Robust-Principal-Component-Analysis-TRPCA" class="textlink">Tensor robust principal component analysis </a></li>       
  <li> <a href="https://github.com/canyilu/tensor-completion-tensor-recovery" class="textlink">Low tubal tensor completion and tensor recovery from Gaussian measurements </a></li>
  <li> <a href="https://github.com/canyilu/LibADMM" class="textlink">A Library of ADMM for Sparse and Low-rank Optimization </a></li>
</ul>

### References
<table class="docutils footnote" frame="void" id="id2" rules="none">
<colgroup><col class="label" /><col /></colgroup>
<tbody valign="top">
<tr><td class="label"><a class="fn-backref" href="#id2">[1]</a></td><td>M. E. Kilmer and C. D. Martin. Factorization strategies for third-order tensors. Linear Algebra and its Applications. 435(3):641–658, 2011.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[2]</a></td><td>C. Lu, J. Feng, Y. Chen, W. Liu, Z. Lin, and S. Yan. Tensor robust principal component analysis with a new tensor nuclear norm. arXiv preprint arXiv:1804.03728, 2018.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[3]</a></td><td>C. Lu, J. Feng, Y. Chen, W. Liu, Z. Lin, and S. Yan. Tensor robust principal component analysis: Exact recovery of corrupted low-rank tensors via convex optimization. In IEEE International Conference on Computer Vision and Pattern Recognition, 2016.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[4]</a></td><td>C. Lu, J. Feng, Z. Lin, and S. Yan. Exact low tubal rank tensor recovery from Gaussian measurements. In International Joint Conference on Artificial Intelligence, 2018.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[5]</a></td><td>C. Lu, J. Feng, S. Yan, Z. Lin. A Unified Alternating Direction Method of Multipliers by Majorization Minimization. IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 40, pp. 527-541, 2018.
</td></tr>
</tbody>
</table>




