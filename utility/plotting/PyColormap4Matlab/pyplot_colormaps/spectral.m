function spectral = spectral()
spectral = [ ... 
       6.196078e-01 3.921569e-03 2.588235e-01 
       6.280661e-01 1.330258e-02 2.608228e-01 
       6.365244e-01 2.268358e-02 2.628220e-01 
       6.449827e-01 3.206459e-02 2.648212e-01 
       6.534410e-01 4.144560e-02 2.668205e-01 
       6.618993e-01 5.082661e-02 2.688197e-01 
       6.703576e-01 6.020761e-02 2.708189e-01 
       6.788158e-01 6.958862e-02 2.728181e-01 
       6.872741e-01 7.896963e-02 2.748174e-01 
       6.957324e-01 8.835063e-02 2.768166e-01 
       7.041907e-01 9.773164e-02 2.788158e-01 
       7.126490e-01 1.071126e-01 2.808151e-01 
       7.211073e-01 1.164937e-01 2.828143e-01 
       7.295656e-01 1.258747e-01 2.848135e-01 
       7.380238e-01 1.352557e-01 2.868128e-01 
       7.464821e-01 1.446367e-01 2.888120e-01 
       7.549404e-01 1.540177e-01 2.908112e-01 
       7.633987e-01 1.633987e-01 2.928105e-01 
       7.718570e-01 1.727797e-01 2.948097e-01 
       7.803153e-01 1.821607e-01 2.968089e-01 
       7.887735e-01 1.915417e-01 2.988082e-01 
       7.972318e-01 2.009227e-01 3.008074e-01 
       8.056901e-01 2.103037e-01 3.028066e-01 
       8.141484e-01 2.196847e-01 3.048058e-01 
       8.226067e-01 2.290657e-01 3.068051e-01 
       8.310650e-01 2.384468e-01 3.088043e-01 
       8.376778e-01 2.467512e-01 3.088812e-01 
       8.424452e-01 2.539792e-01 3.070358e-01 
       8.472126e-01 2.612072e-01 3.051903e-01 
       8.519800e-01 2.684352e-01 3.033449e-01 
       8.567474e-01 2.756632e-01 3.014994e-01 
       8.615148e-01 2.828912e-01 2.996540e-01 
       8.662822e-01 2.901192e-01 2.978085e-01 
       8.710496e-01 2.973472e-01 2.959631e-01 
       8.758170e-01 3.045752e-01 2.941176e-01 
       8.805844e-01 3.118032e-01 2.922722e-01 
       8.853518e-01 3.190311e-01 2.904268e-01 
       8.901192e-01 3.262591e-01 2.885813e-01 
       8.948866e-01 3.334871e-01 2.867359e-01 
       8.996540e-01 3.407151e-01 2.848904e-01 
       9.044214e-01 3.479431e-01 2.830450e-01 
       9.091888e-01 3.551711e-01 2.811995e-01 
       9.139562e-01 3.623991e-01 2.793541e-01 
       9.187236e-01 3.696271e-01 2.775087e-01 
       9.234910e-01 3.768551e-01 2.756632e-01 
       9.282584e-01 3.840830e-01 2.738178e-01 
       9.330258e-01 3.913110e-01 2.719723e-01 
       9.377932e-01 3.985390e-01 2.701269e-01 
       9.425606e-01 4.057670e-01 2.682814e-01 
       9.473280e-01 4.129950e-01 2.664360e-01 
       9.520953e-01 4.202230e-01 2.645905e-01 
       9.568627e-01 4.274510e-01 2.627451e-01 
       9.582468e-01 4.374471e-01 2.673587e-01 
       9.596309e-01 4.474433e-01 2.719723e-01 
       9.610150e-01 4.574394e-01 2.765859e-01 
       9.623991e-01 4.674356e-01 2.811995e-01 
       9.637832e-01 4.774318e-01 2.858131e-01 
       9.651672e-01 4.874279e-01 2.904268e-01 
       9.665513e-01 4.974241e-01 2.950404e-01 
       9.679354e-01 5.074202e-01 2.996540e-01 
       9.693195e-01 5.174164e-01 3.042676e-01 
       9.707036e-01 5.274125e-01 3.088812e-01 
       9.720877e-01 5.374087e-01 3.134948e-01 
       9.734717e-01 5.474048e-01 3.181084e-01 
       9.748558e-01 5.574010e-01 3.227220e-01 
       9.762399e-01 5.673972e-01 3.273356e-01 
       9.776240e-01 5.773933e-01 3.319493e-01 
       9.790081e-01 5.873895e-01 3.365629e-01 
       9.803922e-01 5.973856e-01 3.411765e-01 
       9.817762e-01 6.073818e-01 3.457901e-01 
       9.831603e-01 6.173779e-01 3.504037e-01 
       9.845444e-01 6.273741e-01 3.550173e-01 
       9.859285e-01 6.373702e-01 3.596309e-01 
       9.873126e-01 6.473664e-01 3.642445e-01 
       9.886967e-01 6.573626e-01 3.688581e-01 
       9.900807e-01 6.673587e-01 3.734717e-01 
       9.914648e-01 6.773549e-01 3.780854e-01 
       9.922338e-01 6.861976e-01 3.836217e-01 
       9.923875e-01 6.938870e-01 3.900807e-01 
       9.925413e-01 7.015763e-01 3.965398e-01 
       9.926951e-01 7.092657e-01 4.029988e-01 
       9.928489e-01 7.169550e-01 4.094579e-01 
       9.930027e-01 7.246444e-01 4.159170e-01 
       9.931565e-01 7.323337e-01 4.223760e-01 
       9.933103e-01 7.400231e-01 4.288351e-01 
       9.934641e-01 7.477124e-01 4.352941e-01 
       9.936178e-01 7.554018e-01 4.417532e-01 
       9.937716e-01 7.630911e-01 4.482122e-01 
       9.939254e-01 7.707805e-01 4.546713e-01 
       9.940792e-01 7.784698e-01 4.611303e-01 
       9.942330e-01 7.861592e-01 4.675894e-01 
       9.943868e-01 7.938485e-01 4.740484e-01 
       9.945406e-01 8.015379e-01 4.805075e-01 
       9.946943e-01 8.092272e-01 4.869666e-01 
       9.948481e-01 8.169166e-01 4.934256e-01 
       9.950019e-01 8.246059e-01 4.998847e-01 
       9.951557e-01 8.322953e-01 5.063437e-01 
       9.953095e-01 8.399846e-01 5.128028e-01 
       9.954633e-01 8.476740e-01 5.192618e-01 
       9.956171e-01 8.553633e-01 5.257209e-01 
       9.957709e-01 8.630527e-01 5.321799e-01 
       9.959246e-01 8.707420e-01 5.386390e-01 
       9.960784e-01 8.784314e-01 5.450980e-01 
       9.962322e-01 8.831988e-01 5.530950e-01 
       9.963860e-01 8.879662e-01 5.610919e-01 
       9.965398e-01 8.927336e-01 5.690888e-01 
       9.966936e-01 8.975010e-01 5.770857e-01 
       9.968474e-01 9.022684e-01 5.850827e-01 
       9.970012e-01 9.070358e-01 5.930796e-01 
       9.971549e-01 9.118032e-01 6.010765e-01 
       9.973087e-01 9.165705e-01 6.090734e-01 
       9.974625e-01 9.213379e-01 6.170704e-01 
       9.976163e-01 9.261053e-01 6.250673e-01 
       9.977701e-01 9.308727e-01 6.330642e-01 
       9.979239e-01 9.356401e-01 6.410611e-01 
       9.980777e-01 9.404075e-01 6.490581e-01 
       9.982314e-01 9.451749e-01 6.570550e-01 
       9.983852e-01 9.499423e-01 6.650519e-01 
       9.985390e-01 9.547097e-01 6.730488e-01 
       9.986928e-01 9.594771e-01 6.810458e-01 
       9.988466e-01 9.642445e-01 6.890427e-01 
       9.990004e-01 9.690119e-01 6.970396e-01 
       9.991542e-01 9.737793e-01 7.050365e-01 
       9.993080e-01 9.785467e-01 7.130334e-01 
       9.994617e-01 9.833141e-01 7.210304e-01 
       9.996155e-01 9.880815e-01 7.290273e-01 
       9.997693e-01 9.928489e-01 7.370242e-01 
       9.999231e-01 9.976163e-01 7.450211e-01 
       9.980777e-01 9.992311e-01 7.460208e-01 
       9.942330e-01 9.976932e-01 7.400231e-01 
       9.903883e-01 9.961553e-01 7.340254e-01 
       9.865436e-01 9.946175e-01 7.280277e-01 
       9.826990e-01 9.930796e-01 7.220300e-01 
       9.788543e-01 9.915417e-01 7.160323e-01 
       9.750096e-01 9.900038e-01 7.100346e-01 
       9.711649e-01 9.884660e-01 7.040369e-01 
       9.673203e-01 9.869281e-01 6.980392e-01 
       9.634756e-01 9.853902e-01 6.920415e-01 
       9.596309e-01 9.838524e-01 6.860438e-01 
       9.557862e-01 9.823145e-01 6.800461e-01 
       9.519416e-01 9.807766e-01 6.740484e-01 
       9.480969e-01 9.792388e-01 6.680507e-01 
       9.442522e-01 9.777009e-01 6.620531e-01 
       9.404075e-01 9.761630e-01 6.560554e-01 
       9.365629e-01 9.746251e-01 6.500577e-01 
       9.327182e-01 9.730873e-01 6.440600e-01 
       9.288735e-01 9.715494e-01 6.380623e-01 
       9.250288e-01 9.700115e-01 6.320646e-01 
       9.211842e-01 9.684737e-01 6.260669e-01 
       9.173395e-01 9.669358e-01 6.200692e-01 
       9.134948e-01 9.653979e-01 6.140715e-01 
       9.096501e-01 9.638601e-01 6.080738e-01 
       9.058055e-01 9.623222e-01 6.020761e-01 
       9.019608e-01 9.607843e-01 5.960784e-01 
       8.928874e-01 9.570934e-01 5.979239e-01 
       8.838139e-01 9.534025e-01 5.997693e-01 
       8.747405e-01 9.497116e-01 6.016148e-01 
       8.656671e-01 9.460208e-01 6.034602e-01 
       8.565936e-01 9.423299e-01 6.053057e-01 
       8.475202e-01 9.386390e-01 6.071511e-01 
       8.384468e-01 9.349481e-01 6.089965e-01 
       8.293733e-01 9.312572e-01 6.108420e-01 
       8.202999e-01 9.275663e-01 6.126874e-01 
       8.112265e-01 9.238754e-01 6.145329e-01 
       8.021530e-01 9.201845e-01 6.163783e-01 
       7.930796e-01 9.164937e-01 6.182238e-01 
       7.840062e-01 9.128028e-01 6.200692e-01 
       7.749327e-01 9.091119e-01 6.219146e-01 
       7.658593e-01 9.054210e-01 6.237601e-01 
       7.567859e-01 9.017301e-01 6.256055e-01 
       7.477124e-01 8.980392e-01 6.274510e-01 
       7.386390e-01 8.943483e-01 6.292964e-01 
       7.295656e-01 8.906574e-01 6.311419e-01 
       7.204921e-01 8.869666e-01 6.329873e-01 
       7.114187e-01 8.832757e-01 6.348328e-01 
       7.023453e-01 8.795848e-01 6.366782e-01 
       6.932718e-01 8.758939e-01 6.385236e-01 
       6.841984e-01 8.722030e-01 6.403691e-01 
       6.751250e-01 8.685121e-01 6.422145e-01 
       6.652826e-01 8.645905e-01 6.432141e-01 
       6.546713e-01 8.604383e-01 6.433679e-01 
       6.440600e-01 8.562860e-01 6.435217e-01 
       6.334487e-01 8.521338e-01 6.436755e-01 
       6.228374e-01 8.479815e-01 6.438293e-01 
       6.122261e-01 8.438293e-01 6.439831e-01 
       6.016148e-01 8.396770e-01 6.441369e-01 
       5.910035e-01 8.355248e-01 6.442907e-01 
       5.803922e-01 8.313725e-01 6.444444e-01 
       5.697809e-01 8.272203e-01 6.445982e-01 
       5.591696e-01 8.230681e-01 6.447520e-01 
       5.485582e-01 8.189158e-01 6.449058e-01 
       5.379469e-01 8.147636e-01 6.450596e-01 
       5.273356e-01 8.106113e-01 6.452134e-01 
       5.167243e-01 8.064591e-01 6.453672e-01 
       5.061130e-01 8.023068e-01 6.455210e-01 
       4.955017e-01 7.981546e-01 6.456747e-01 
       4.848904e-01 7.940023e-01 6.458285e-01 
       4.742791e-01 7.898501e-01 6.459823e-01 
       4.636678e-01 7.856978e-01 6.461361e-01 
       4.530565e-01 7.815456e-01 6.462899e-01 
       4.424452e-01 7.773933e-01 6.464437e-01 
       4.318339e-01 7.732411e-01 6.465975e-01 
       4.212226e-01 7.690888e-01 6.467512e-01 
       4.106113e-01 7.649366e-01 6.469050e-01 
       4.000000e-01 7.607843e-01 6.470588e-01 
       3.920031e-01 7.518647e-01 6.507497e-01 
       3.840062e-01 7.429450e-01 6.544406e-01 
       3.760092e-01 7.340254e-01 6.581315e-01 
       3.680123e-01 7.251057e-01 6.618224e-01 
       3.600154e-01 7.161861e-01 6.655133e-01 
       3.520185e-01 7.072664e-01 6.692042e-01 
       3.440215e-01 6.983468e-01 6.728950e-01 
       3.360246e-01 6.894271e-01 6.765859e-01 
       3.280277e-01 6.805075e-01 6.802768e-01 
       3.200308e-01 6.715879e-01 6.839677e-01 
       3.120338e-01 6.626682e-01 6.876586e-01 
       3.040369e-01 6.537486e-01 6.913495e-01 
       2.960400e-01 6.448289e-01 6.950404e-01 
       2.880431e-01 6.359093e-01 6.987313e-01 
       2.800461e-01 6.269896e-01 7.024221e-01 
       2.720492e-01 6.180700e-01 7.061130e-01 
       2.640523e-01 6.091503e-01 7.098039e-01 
       2.560554e-01 6.002307e-01 7.134948e-01 
       2.480584e-01 5.913110e-01 7.171857e-01 
       2.400615e-01 5.823914e-01 7.208766e-01 
       2.320646e-01 5.734717e-01 7.245675e-01 
       2.240677e-01 5.645521e-01 7.282584e-01 
       2.160707e-01 5.556324e-01 7.319493e-01 
       2.080738e-01 5.467128e-01 7.356401e-01 
       2.000769e-01 5.377932e-01 7.393310e-01 
       1.994617e-01 5.289504e-01 7.391003e-01 
       2.062284e-01 5.201845e-01 7.349481e-01 
       2.129950e-01 5.114187e-01 7.307958e-01 
       2.197616e-01 5.026528e-01 7.266436e-01 
       2.265283e-01 4.938870e-01 7.224913e-01 
       2.332949e-01 4.851211e-01 7.183391e-01 
       2.400615e-01 4.763552e-01 7.141869e-01 
       2.468281e-01 4.675894e-01 7.100346e-01 
       2.535948e-01 4.588235e-01 7.058824e-01 
       2.603614e-01 4.500577e-01 7.017301e-01 
       2.671280e-01 4.412918e-01 6.975779e-01 
       2.738947e-01 4.325260e-01 6.934256e-01 
       2.806613e-01 4.237601e-01 6.892734e-01 
       2.874279e-01 4.149942e-01 6.851211e-01 
       2.941945e-01 4.062284e-01 6.809689e-01 
       3.009612e-01 3.974625e-01 6.768166e-01 
       3.077278e-01 3.886967e-01 6.726644e-01 
       3.144944e-01 3.799308e-01 6.685121e-01 
       3.212611e-01 3.711649e-01 6.643599e-01 
       3.280277e-01 3.623991e-01 6.602076e-01 
       3.347943e-01 3.536332e-01 6.560554e-01 
       3.415609e-01 3.448674e-01 6.519031e-01 
       3.483276e-01 3.361015e-01 6.477509e-01 
       3.550942e-01 3.273356e-01 6.435986e-01 
       3.618608e-01 3.185698e-01 6.394464e-01 
       3.686275e-01 3.098039e-01 6.352941e-01 ] ; 
