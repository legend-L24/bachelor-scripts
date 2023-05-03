Theory
======

To predict charges with QEq, we multiply the inverse hardness matrix :math:`\bar{\mathbf{A}}` with the vector of electronegativities :math:`\mathbf{c}`:

.. math::

    \mathbf{q} = - \bar{\mathbf{A}} \mathbf{c}



We now substitute the vector of electronegativities with a KRR expression:

.. math::

    \mathbf{c} = \mathbf{K} \alpha


so that

.. math::

    \mathbf{q} = - \bar{\mathbf{A}} \mathbf{K} \alpha


Finally, we can calculate the dipole vector of a molecule using the :math:`3 \times N` matrix :math:`\mathbf{R}` of center-of-mass shifted atomic coordinates:

.. math::

    \mathbf{\mu} = - \mathbf{R} (\bar{\mathbf{A}} \mathbf{K} \alpha) 

We now want to determine the optimal regression weights :math:`\alpha`. To this end we set up a regularized least-squares problem:

.. math::

    L = ||(\mathbf{\mu}-\mathbf{\mu_{ref}})||^2 + \lambda \alpha^T \mathbf{K} \alpha = || - \mathbf{R} (\bar{\mathbf{A}} \mathbf{K} \alpha)  -\mathbf{\mu_{ref}}||^2 + \lambda \alpha^T \mathbf{K} \alpha


To minimize this loss function, we set :math:`\frac{dL}{d\alpha}=0` and solve for :math:`\alpha`, to obtain: 

.. math::
    
    \alpha = - (\bar{\mathbf{A}}^T \mathbf{R}^T \mathbf{R} \bar{\mathbf{A}} \mathbf{K} + \lambda \mathbb{1})^{-1} \bar{\mathbf{A}}^T \mathbf{R}^T \mathbf{\mu_{ref}}


The above equations are formulated for a single QEq problem (i.e. a single molecule or system). In practice we train on multiple systems at once. This can still be achieved in a single linear algebra equation by using blocked matrices for :math:`\bar{\mathbf{A}}` and :math:`\mathbf{R}`, and concatenating the dipole vector elements of all training systems into a single vector of dimension :math:`3N`. 

Apart from the electronegativities, QEq models depend on several other parameters, which are set via simple heuristics: To model the charge density, Gaussian charge distributions are used. The standard deviation of these distributions is assumed to be proportional to the original QEq radii of the elements. This scaling can be specified via the parameter :code:`scale_atsize`, where the default value of 1.0 works well for the cases we tested. The elemental hardness is directly derived from the on-site term of the electrostatic energy of the Gaussians, so there's no additional parameter used here. The only other hyperparameter of the model is the regularization stength :math:`\lambda`. This means that two global parameter (:code:`scale_atsize` and :code:`lambda_reg`) are enough to completely specify a kQEq model. Of course, it would be straightforward to use separate atomic radii or independent electronic hardness parameters, if more flexibility is required (this is also implemented).

For a detailed derivation of the kQEq equations please refer to !paper!
