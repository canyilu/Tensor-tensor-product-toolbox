## Tensor-Tensor Product Toolbox

### Introduction

The tensor-tensor product (t-product) <a class="footnote-reference" href="#id2" id="id1">[1]</a> is a natural generalization of matrix multiplication. Based on t-product, many operations on matrix can be extended to tensor cases, including tensor SVD (see an illustration in the figure below), tensor spectral norm, tensor nuclear norm <a class="footnote-reference" href="#id2" id="id1">[2]</a> and many others. The linear algebraic structure of tensors are similar to the matrix cases. We develop a Matlab toolbox to implement several basic operations on tensors based on t-product.

<p align="center"> 
<img src="https://github.com/canyilu/tproduct/blob/master/tsvd.JPG">
</p>

### List of functions

The table below gives the list of functions implemented in our toolbox. The detailed definitions of these tensor concepts, operations and tensor factorizations are given at <a href="../tproduct/blob/master/manual.pdf" class="textlink" target="_blank">https://github.com/canyilu/tproduct/blob/master/manual.pdf</a>. 

<p align="center"> 
<img src="https://github.com/canyilu/tproduct/blob/master/tab_funs.JPG">
</p>

Note that we only focus on 3 way tensor in this toolbox. We will develop the same functions for p-way tensor in the near future. We will also provide the python version soon.

Simply run the following file to test all the above functions:
```matlab
test.m
```

### Citing

<p>In citing this toolbox in your papers, please use the following reference:</p>

<blockquote></blockquote>
<div><p>Canyi Lu, Tensor-Tensor Product Toolbox. Carnegie Mellon University. June, 2018.
<tt class="docutils literal"><span class="pre">https://github.com/canyilu/tproduct</span></tt>.</p>
</div>

<p>The corresponding BiBTeX citation are given below:</p>
<div class="highlight-none"><div class="highlight"><pre>
@incollection{tproduct2018lu,
  author    = {Lu, Canyi},
  title     = {Tensor-Tensor Product Toolbox},
  publisher = {Carnegie Mellon University},
  month     = {June},
  year      = {2018},
  note      = {\url{https://github.com/canyilu/tproduct}}
}
</pre></div>
  
  
### References
<table class="docutils footnote" frame="void" id="id2" rules="none">
<colgroup><col class="label" /><col /></colgroup>
<tbody valign="top">
<tr><td class="label"><a class="fn-backref" href="#id2">[1]</a></td><td>M. E. Kilmer and C. D. Martin. Factorization strategies for third-order tensors. Linear Algebra and its Applications, 435(3):641–658, 2011.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[2]</a></td><td>C. Lu, J. Feng, Y. Chen, W. Liu, Z. Lin, and S. Yan. Tensor robust principal component analysis with a new tensor nuclear norm. arXiv preprint arXiv:1804.03728, 2018.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[3]</a></td><td>C. Lu, J. Feng, Y. Chen, W. Liu, Z. Lin, and S. Yan. Tensor robust principal component analysis: Exact recovery of corrupted low-rank tensors via convex optimization. In IEEE International Conference on Computer Vision and Pattern Recognition, 2016.</td></tr>
<tr><td class="label"><a class="fn-backref" href="#id2">[4]</a></td><td>C. Lu, J. Feng, Z. Lin, and S. Yan. Exact low tubal rank tensor recovery from gaussian measurements. In International Joint Conference on Artificial Intelligence, 2018.
</td></tr>
</tbody>
</table>




