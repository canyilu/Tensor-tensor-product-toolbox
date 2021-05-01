## Tensor-Tensor Product Toolbox 2.0 (updated in April, 2021)

### 1. T-product Toolbox 1.0

The tensor-tensor product (t-product) <a class="footnote-reference" href="#id2" id="id1">[1]</a> is a natural generalization of matrix multiplication. Based on t-product, many operations on matrix can be extended to tensor cases, including tensor SVD (see an illustration in the figure below), tensor spectral norm, tensor nuclear norm <a class="footnote-reference" href="#id2" id="id1">[2]</a> and many others. 
<p align="center"> 
<img src="https://github.com/canyilu/tproduct/blob/master/doc/figure_tsvd.JPG">
</p>

The linear algebraic structure of tensors are similar to the matrix cases. We have tensor-tensor product, tensor SVD, tensor inverse and some other reated concepts extended from matrices. The detailed definitions of these tensor concepts, operations and tensor factorizations are given at <a href="../publications/2018-software-tproduct.pdf" class="textlink" target="_blank">https://canyilu.github.io/publications/2018-software-tproduct.pdf</a>. We develop a Matlab toolbox to implement several basic operations on tensors based on t-product. See a list of implemented functions in t-product toolbox 1.0 below.
<p align="center"> 
<img src="https://github.com/canyilu/tproduct/blob/master/doc/figure_functions_tproduct_1.0.JPG">
</p>

### 2. T-product Toolbox 2.0

The original t-product <a class="footnote-reference" href="#id2" id="id1">[1]</a> uses the discrete Fourier transform and uses the fast Fourier transform (FFT) for efficient computing. It is further generlaized to the t-product under arbitrary invertible linear transform in <a class="footnote-reference" href="#id2" id="id1">[2]</a>. Thus, all the concepts (e.g., tsvd, tensor inverse) of t-product under FFT can be generalized to t-product under general linear transforms. If the linear transform satisfies <a href="https://www.codecogs.com/eqnedit.php?latex=L^\top&space;L&space;=&space;LL^\top=\ell&space;I" target="_blank"><img src="https://latex.codecogs.com/gif.latex?L^\top&space;L&space;=&space;LL^\top=\ell&space;I" title="L^\top L = LL^\top=\ell I" /></a> where <a href="https://www.codecogs.com/eqnedit.php?latex=\ell>0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell>0" title="\ell>0" /></a>, then we can define a more general tensor nuclear norm induced by the t-product under this linear transform. Then we develop a more general Matlab toolbox to implement t-product under general linear transform. See a list of implemented functions in t-product toolbox 2.0 below.
<p align="center"> 
<img src="https://github.com/canyilu/tproduct/blob/master/doc/figure_functions_tproduct_2.0.JPG">
</p>

For the definitions of t-product and related concepts under linear transform, please refer to <a class="footnote-reference" href="#id2" id="id1">[2]</a> and our works <a class="footnote-reference" href="#id2" id="id1">[6,7]</a>. We will provide a document to give the details in the future.

### 3. Examples
Simply run the following routine to test all the above functions:
```matlab
test.m
```

### 4. Citation

<p>In citing this toolbox in your papers, please use the following references:</p>

<div class="highlight-none"><div class="highlight"><pre>
C. Lu. Tensor-Tensor Product Toolbox. Carnegie Mellon University, June 2018. https://github.com/canyilu/tproduct.
C. Lu, J. Feng, Y. Chen, W. Liu, Z. Lin, and S. Yan. Tensor robust principal component analysis with a new tensor
nuclear norm. IEEE Transactions on Pattern Analysis and Machine Intelligence, 2019.
C. Lu, X. Peng, and Y. Wei. Low-Rank Tensor Completion With a New Tensor Nuclear Norm Induced by Invertible Linear Transforms. IEEE International Conference on Computer Vision and Pattern Recognition (CVPR), 2019
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
  journal      = {IEEE Transactions on Pattern Analysis and Machine Intelligence},
  year         = {2019}
}
@inproceedings{lu2019tensor,
  author       = {Lu, Canyi and Peng, Xi and Wei, Yunchao},
  title        = {Low-Rank Tensor Completion With a New Tensor Nuclear Norm Induced by Invertible Linear Transforms},
  journal      = {CVPR},
  year         = {2019}
}
</pre></div>


### 5. Version History
- Version 1.0 was released on June, 2018. It implements the functions of t-product and related concepts under fast Fourier transform.
- Version 2.0 was released on April, 2021. It implements the functions of t-product and related concepts under general invertible linear transform. The fast Fourier transform is the default transform.
  + Most functions are direct generalization from the fast Fourier transform to general linear transform, e.g., ```tprod```, ```tran```, ```teye```, ```tinv```, ```tsvd```, ```tubalrank```, ```tsn```, ```tnn```, ```prox_tnn``` and ```tqr```.
  + Some functions are new (not included in Version 1.0), e.g., ```basis_column```, ```basis_tube``` and ```unit_eijk```.
  + Some functions in Version 1.0 are updated, e.g., the setting of parameter tol in ```tubalrank``` and ```tsvd``` is updated, and ```tprod```, ```tsn```, ```tinv``` and ```tqr``` are updated.



### 6. Related Toolboxes
The t-product toolbox has been applied in our works for tensor roubst PCA <a class="footnote-reference" href="#id2" id="id1">[3,4]</a>, low-rank tensor completion and low-rank tensor recovery from Gaussian measurements <a class="footnote-reference" href="#id2" id="id1">[5]</a>. The t-product under linear transform has also been applied in tensor completion <a class="footnote-reference" href="#id2" id="id1">[6]</a> and tensor robust PCA <a class="footnote-reference" href="#id2" id="id1">[7]</a>. Some more models are included in LibADMM toolbox <a class="footnote-reference" href="#id2" id="id1">[8]</a>.
<ul>
  <li> <a href="https://github.com/canyilu/Tensor-Robust-Principal-Component-Analysis-TRPCA" class="textlink">Tensor robust principal component analysis </a></li>       
  <li> <a href="https://github.com/canyilu/tensor-completion-tensor-recovery" class="textlink">Low tubal tensor completion and tensor recovery from Gaussian measurements </a></li>
  <li> <a href="https://github.com/canyilu/Tensor-robust-PCA-and-tensor-completion-under-linear-transform" class="textlink">Tensor robust PCA and tensor completion based on tensor nuclear norm under linear transform</a></li>
  <li> <a href="https://github.com/canyilu/LibADMM" class="textlink">A Library of ADMM for Sparse and Low-rank Optimization </a></li>
</ul>

References
<table class="docutils footnote" frame="void" id="id2" rules="none">
<colgroup><col class="label" /><col /></colgroup>
<tbody valign="top">
<tr><td class="label"><a class="fn-backref" href="#id2">[1]</a></td><td>M. E. Kilmer and C. D. Martin. Factorization strategies for third-order tensors. Linear Algebra and its Applications. 435(3):641–658, 2011.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[2]</a></td><td>M. E. Kilmer and S. Aeron. Tensor-Tensor Products with Invertible Linear Transforms. Linear Algebra and its Applications. 2015.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[3]</a></td><td>C. Lu, J. Feng, Y. Chen, W. Liu, Z. Lin, and S. Yan. Tensor robust principal component analysis with a new tensor nuclear norm. IEEE Transactions on Pattern Analysis and Machine Intelligence, 2019.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[4]</a></td><td>C. Lu, J. Feng, Y. Chen, W. Liu, Z. Lin, and S. Yan. Tensor robust principal component analysis: Exact recovery of corrupted low-rank tensors via convex optimization. In IEEE International Conference on Computer Vision and Pattern Recognition, 2016.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[5]</a></td><td>C. Lu, J. Feng, Z. Lin, and S. Yan. Exact low tubal rank tensor recovery from Gaussian measurements. In International Joint Conference on Artificial Intelligence, 2018.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[6]</a></td><td>C. Lu, X. Peng, and Y. Wei. Low-Rank Tensor Completion With a New Tensor Nuclear Norm Induced by Invertible Linear Transforms. IEEE International Conference on Computer Vision and Pattern Recognition (CVPR), 2019.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[7]</a></td><td>C. Lu. Exact Recovery of Tensor Robust Principal Component Analysis under Linear Transforms. arXiv preprint arXiv:1907.08288. 2019.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[8]</a></td><td>C. Lu, J. Feng, S. Yan, Z. Lin. A Unified Alternating Direction Method of Multipliers by Majorization Minimization. IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 40, pp. 527-541, 2018.
</td></tr>
</tbody>
</table>




