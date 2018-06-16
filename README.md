# Introduction

Tensors are higher-order extensions of matrices. In recent work <a class="footnote-reference" href="#id2" id="id1">[1]</a>, the authors introduced the notion of the t-product, a generalization of matrix multiplication for tensors of order three. The multiplication is based on a convolution-like operation, which 	can be implemented efficiently using the Fast Fourier Transform (FFT). Based on t-product, there has a similar linear algebraic structure of tensors to matrices. For example, there has the tensor SVD (t-SVD) which is computable. By using some properties of FFT, we have a more efficient way for computing t-product and t-SVD in  <a class="footnote-reference" href="#id2" id="id1">[2]</a>. We develop a Matlab toolbox to implement several basic operations on tensors based on t-product.
![Alt text](https://github.com/canyilu/tproduct/blob/master/tsvd.JPG)


The table below gives the list of functions implemented in our toolbox. The detailed definitions of these tensor concepts, operations and tensor factorizations are given at <a href="../tproduct/manual.pdf" class="textlink" target="_blank">https://github.com/canyilu/tproduct/manual.pdf</a>. Note that we only focus on 3 way tensor in this toolbox. We will develop the same functions for p-way tensor in the near future. We will also provide the python verion soon.
![Alt text](https://github.com/canyilu/tproduct/blob/master/tab_tprod_funlist.JPG)

# Citation

<p>In citing this toolbox in your papers, please use the following reference:</p>

<blockquote>
<div><p>Canyi Lu, Tensor-Tensor Product Toolbox. Carnegie Mellon University. June, 2018.
<tt class="docutils literal"><span class="pre">https://github.com/canyilu/tproduct</span></tt>.</p>
</div></blockquote>

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
  
  
# Reference
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




