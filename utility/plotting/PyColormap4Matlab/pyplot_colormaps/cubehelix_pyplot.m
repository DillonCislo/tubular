function cubehelix = cubehelix_pyplot()
cubehelix = [ ... 
       0 0 0 
       6.716295e-03 2.118574e-03 5.970233e-03 
       1.325242e-02 4.287499e-03 1.216179e-02 
       1.959920e-02 6.513601e-03 1.856304e-02 
       2.574810e-02 8.803482e-03 2.516187e-02 
       3.169123e-02 1.116350e-02 3.194570e-02 
       3.742134e-02 1.359977e-02 3.890155e-02 
       4.293183e-02 1.611812e-02 4.601605e-02 
       4.821682e-02 1.872413e-02 5.327545e-02 
       5.327106e-02 2.142305e-02 6.066574e-02 
       5.809005e-02 2.421988e-02 6.817258e-02 
       6.266993e-02 2.711926e-02 7.578143e-02 
       6.700758e-02 3.012554e-02 8.347751e-02 
       7.110058e-02 3.324274e-02 9.124592e-02 
       7.494720e-02 3.647454e-02 9.907158e-02 
       7.854643e-02 3.982428e-02 1.069394e-01 
       8.189796e-02 4.329495e-02 1.148340e-01 
       8.500218e-02 4.688919e-02 1.227404e-01 
       8.786018e-02 5.060929e-02 1.306433e-01 
       9.047371e-02 5.445716e-02 1.385275e-01 
       9.284524e-02 5.843437e-02 1.463780e-01 
       9.497788e-02 6.254212e-02 1.541799e-01 
       9.687542e-02 6.678122e-02 1.619184e-01 
       9.854228e-02 7.115216e-02 1.695789e-01 
       9.998353e-02 7.565501e-02 1.771472e-01 
       1.012048e-01 8.028952e-02 1.846092e-01 
       1.022125e-01 8.505507e-02 1.919510e-01 
       1.030133e-01 8.995065e-02 1.991594e-01 
       1.036148e-01 9.497494e-02 2.062211e-01 
       1.040249e-01 1.001262e-01 2.131236e-01 
       1.042521e-01 1.054025e-01 2.198546e-01 
       1.043053e-01 1.108014e-01 2.264022e-01 
       1.041942e-01 1.163202e-01 2.327552e-01 
       1.039285e-01 1.219558e-01 2.389027e-01 
       1.035188e-01 1.277050e-01 2.448344e-01 
       1.029756e-01 1.335640e-01 2.505406e-01 
       1.023103e-01 1.395290e-01 2.560120e-01 
       1.015341e-01 1.455956e-01 2.612401e-01 
       1.006591e-01 1.517594e-01 2.662170e-01 
       9.969721e-02 1.580155e-01 2.709351e-01 
       9.866087e-02 1.643590e-01 2.753878e-01 
       9.756267e-02 1.707846e-01 2.795691e-01 
       9.641541e-02 1.772867e-01 2.834736e-01 
       9.523206e-02 1.838597e-01 2.870966e-01 
       9.402576e-02 1.904975e-01 2.904340e-01 
       9.280972e-02 1.971942e-01 2.934826e-01 
       9.159727e-02 2.039434e-01 2.962397e-01 
       9.040175e-02 2.107386e-01 2.987034e-01 
       8.923656e-02 2.175734e-01 3.008727e-01 
       8.811505e-02 2.244409e-01 3.027469e-01 
       8.705056e-02 2.313344e-01 3.043265e-01 
       8.605634e-02 2.382469e-01 3.056124e-01 
       8.514552e-02 2.451715e-01 3.066062e-01 
       8.433111e-02 2.521010e-01 3.073104e-01 
       8.362597e-02 2.590284e-01 3.077282e-01 
       8.304274e-02 2.659465e-01 3.078632e-01 
       8.259386e-02 2.728481e-01 3.077201e-01 
       8.229149e-02 2.797261e-01 3.073040e-01 
       8.214755e-02 2.865734e-01 3.066208e-01 
       8.217361e-02 2.933829e-01 3.056770e-01 
       8.238094e-02 3.001474e-01 3.044797e-01 
       8.278044e-02 3.068600e-01 3.030367e-01 
       8.338263e-02 3.135138e-01 3.013563e-01 
       8.419763e-02 3.201020e-01 2.994476e-01 
       8.523512e-02 3.266178e-01 2.973201e-01 
       8.650432e-02 3.330547e-01 2.949839e-01 
       8.801402e-02 3.394063e-01 2.924495e-01 
       8.977246e-02 3.456663e-01 2.897281e-01 
       9.178742e-02 3.518285e-01 2.868312e-01 
       9.406612e-02 3.578871e-01 2.837710e-01 
       9.661525e-02 3.638364e-01 2.805598e-01 
       9.944094e-02 3.696707e-01 2.772105e-01 
       1.025487e-01 3.753849e-01 2.737365e-01 
       1.059436e-01 3.809739e-01 2.701511e-01 
       1.096299e-01 3.864328e-01 2.664684e-01 
       1.136115e-01 3.917572e-01 2.627025e-01 
       1.178913e-01 3.969426e-01 2.588678e-01 
       1.224721e-01 4.019851e-01 2.549790e-01 
       1.273556e-01 4.068811e-01 2.510510e-01 
       1.325432e-01 4.116269e-01 2.470986e-01 
       1.380353e-01 4.162196e-01 2.431371e-01 
       1.438321e-01 4.206563e-01 2.391816e-01 
       1.499328e-01 4.249345e-01 2.352474e-01 
       1.563362e-01 4.290519e-01 2.313498e-01 
       1.630402e-01 4.330068e-01 2.275041e-01 
       1.700423e-01 4.367976e-01 2.237256e-01 
       1.773394e-01 4.404231e-01 2.200294e-01 
       1.849274e-01 4.438824e-01 2.164306e-01 
       1.928021e-01 4.471751e-01 2.129442e-01 
       2.009584e-01 4.503009e-01 2.095849e-01 
       2.093905e-01 4.532599e-01 2.063674e-01 
       2.180923e-01 4.560528e-01 2.033060e-01 
       2.270568e-01 4.586803e-01 2.004147e-01 
       2.362769e-01 4.611437e-01 1.977073e-01 
       2.457444e-01 4.634444e-01 1.951974e-01 
       2.554510e-01 4.655842e-01 1.928981e-01 
       2.653876e-01 4.675655e-01 1.908221e-01 
       2.755449e-01 4.693906e-01 1.889817e-01 
       2.859129e-01 4.710625e-01 1.873889e-01 
       2.964811e-01 4.725842e-01 1.860552e-01 
       3.072388e-01 4.739592e-01 1.849915e-01 
       3.181748e-01 4.751913e-01 1.842083e-01 
       3.292773e-01 4.762846e-01 1.837156e-01 
       3.405344e-01 4.772433e-01 1.835227e-01 
       3.519339e-01 4.780721e-01 1.836385e-01 
       3.634629e-01 4.787759e-01 1.840713e-01 
       3.751087e-01 4.793599e-01 1.848286e-01 
       3.868579e-01 4.798294e-01 1.859176e-01 
       3.986973e-01 4.801901e-01 1.873446e-01 
       4.106130e-01 4.804478e-01 1.891154e-01 
       4.225914e-01 4.806087e-01 1.912351e-01 
       4.346184e-01 4.806789e-01 1.937080e-01 
       4.466800e-01 4.806650e-01 1.965379e-01 
       4.587620e-01 4.805737e-01 1.997279e-01 
       4.708500e-01 4.804116e-01 2.032802e-01 
       4.829300e-01 4.801859e-01 2.071966e-01 
       4.949874e-01 4.799035e-01 2.114779e-01 
       5.070081e-01 4.795717e-01 2.161244e-01 
       5.189779e-01 4.791978e-01 2.211357e-01 
       5.308826e-01 4.787892e-01 2.265105e-01 
       5.427083e-01 4.783534e-01 2.322469e-01 
       5.544409e-01 4.778980e-01 2.383424e-01 
       5.660669e-01 4.774304e-01 2.447937e-01 
       5.775726e-01 4.769584e-01 2.515967e-01 
       5.889449e-01 4.764896e-01 2.587468e-01 
       6.001706e-01 4.760316e-01 2.662388e-01 
       6.112369e-01 4.755919e-01 2.740665e-01 
       6.221315e-01 4.751783e-01 2.822234e-01 
       6.328422e-01 4.747981e-01 2.907021e-01 
       6.433573e-01 4.744589e-01 2.994948e-01 
       6.536652e-01 4.741681e-01 3.085929e-01 
       6.637551e-01 4.739329e-01 3.179873e-01 
       6.736163e-01 4.737605e-01 3.276684e-01 
       6.832387e-01 4.736580e-01 3.376260e-01 
       6.926126e-01 4.736323e-01 3.478493e-01 
       7.017290e-01 4.736901e-01 3.583270e-01 
       7.105790e-01 4.738381e-01 3.690474e-01 
       7.191546e-01 4.740827e-01 3.799984e-01 
       7.274482e-01 4.744300e-01 3.911672e-01 
       7.354527e-01 4.748862e-01 4.025409e-01 
       7.431615e-01 4.754570e-01 4.141061e-01 
       7.505688e-01 4.761479e-01 4.258489e-01 
       7.576694e-01 4.769644e-01 4.377552e-01 
       7.644584e-01 4.779115e-01 4.498107e-01 
       7.709317e-01 4.789939e-01 4.620007e-01 
       7.770860e-01 4.802164e-01 4.743104e-01 
       7.829183e-01 4.815830e-01 4.867245e-01 
       7.884265e-01 4.830979e-01 4.992279e-01 
       7.936089e-01 4.847646e-01 5.118051e-01 
       7.984646e-01 4.865865e-01 5.244405e-01 
       8.029934e-01 4.885668e-01 5.371187e-01 
       8.071956e-01 4.907081e-01 5.498239e-01 
       8.110720e-01 4.930129e-01 5.625404e-01 
       8.146245e-01 4.954832e-01 5.752526e-01 
       8.178553e-01 4.981208e-01 5.879449e-01 
       8.207671e-01 5.009272e-01 6.006018e-01 
       8.233636e-01 5.039033e-01 6.132080e-01 
       8.256489e-01 5.070501e-01 6.257480e-01 
       8.276277e-01 5.103678e-01 6.382069e-01 
       8.293053e-01 5.138566e-01 6.505698e-01 
       8.306876e-01 5.175161e-01 6.628221e-01 
       8.317811e-01 5.213458e-01 6.749494e-01 
       8.325929e-01 5.253447e-01 6.869377e-01 
       8.331304e-01 5.295114e-01 6.987732e-01 
       8.334019e-01 5.338445e-01 7.104426e-01 
       8.334159e-01 5.383418e-01 7.219330e-01 
       8.331815e-01 5.430012e-01 7.332317e-01 
       8.327083e-01 5.478200e-01 7.443266e-01 
       8.320062e-01 5.527954e-01 7.552060e-01 
       8.310857e-01 5.579240e-01 7.658589e-01 
       8.299577e-01 5.632024e-01 7.762744e-01 
       8.286333e-01 5.686268e-01 7.864426e-01 
       8.271241e-01 5.741930e-01 7.963538e-01 
       8.254421e-01 5.798968e-01 8.059990e-01 
       8.235994e-01 5.857333e-01 8.153699e-01 
       8.216085e-01 5.916979e-01 8.244585e-01 
       8.194821e-01 5.977853e-01 8.332578e-01 
       8.172332e-01 6.039901e-01 8.417613e-01 
       8.148750e-01 6.103068e-01 8.499629e-01 
       8.124208e-01 6.167296e-01 8.578576e-01 
       8.098839e-01 6.232524e-01 8.654408e-01 
       8.072781e-01 6.298691e-01 8.727086e-01 
       8.046168e-01 6.365734e-01 8.796578e-01 
       8.019139e-01 6.433587e-01 8.862861e-01 
       7.991830e-01 6.502183e-01 8.925917e-01 
       7.964379e-01 6.571457e-01 8.985735e-01 
       7.936921e-01 6.641337e-01 9.042312e-01 
       7.909592e-01 6.711755e-01 9.095651e-01 
       7.882529e-01 6.782641e-01 9.145763e-01 
       7.855863e-01 6.853922e-01 9.192667e-01 
       7.829728e-01 6.925527e-01 9.236386e-01 
       7.804252e-01 6.997385e-01 9.276953e-01 
       7.779565e-01 7.069422e-01 9.314406e-01 
       7.755791e-01 7.141567e-01 9.348791e-01 
       7.733054e-01 7.213747e-01 9.380159e-01 
       7.711472e-01 7.285891e-01 9.408569e-01 
       7.691163e-01 7.357928e-01 9.434086e-01 
       7.672239e-01 7.429786e-01 9.456782e-01 
       7.654809e-01 7.501396e-01 9.476733e-01 
       7.638979e-01 7.572689e-01 9.494023e-01 
       7.624848e-01 7.643597e-01 9.508741e-01 
       7.612515e-01 7.714053e-01 9.520982e-01 
       7.602069e-01 7.783993e-01 9.530844e-01 
       7.593598e-01 7.853352e-01 9.538434e-01 
       7.587183e-01 7.922069e-01 9.543861e-01 
       7.582901e-01 7.990084e-01 9.547239e-01 
       7.580822e-01 8.057338e-01 9.548688e-01 
       7.581011e-01 8.123776e-01 9.548329e-01 
       7.583527e-01 8.189344e-01 9.546289e-01 
       7.588425e-01 8.253990e-01 9.542699e-01 
       7.595750e-01 8.317666e-01 9.537692e-01 
       7.605545e-01 8.380325e-01 9.531403e-01 
       7.617844e-01 8.441925e-01 9.523972e-01 
       7.632674e-01 8.502423e-01 9.515540e-01 
       7.650060e-01 8.561782e-01 9.506249e-01 
       7.670014e-01 8.619968e-01 9.496243e-01 
       7.692548e-01 8.676949e-01 9.485670e-01 
       7.717662e-01 8.732697e-01 9.474674e-01 
       7.745352e-01 8.787185e-01 9.463405e-01 
       7.775608e-01 8.840393e-01 9.452008e-01 
       7.808412e-01 8.892301e-01 9.440632e-01 
       7.843740e-01 8.942894e-01 9.429423e-01 
       7.881561e-01 8.992162e-01 9.418527e-01 
       7.921838e-01 9.040095e-01 9.408089e-01 
       7.964528e-01 9.086690e-01 9.398254e-01 
       8.009580e-01 9.131944e-01 9.389161e-01 
       8.056940e-01 9.175861e-01 9.380950e-01 
       8.106543e-01 9.218447e-01 9.373759e-01 
       8.158323e-01 9.259711e-01 9.367721e-01 
       8.212206e-01 9.299667e-01 9.362966e-01 
       8.268111e-01 9.338331e-01 9.359622e-01 
       8.325955e-01 9.375724e-01 9.357812e-01 
       8.385645e-01 9.411869e-01 9.357656e-01 
       8.447088e-01 9.446794e-01 9.359266e-01 
       8.510181e-01 9.480530e-01 9.362755e-01 
       8.574821e-01 9.513110e-01 9.368225e-01 
       8.640897e-01 9.544571e-01 9.375776e-01 
       8.708296e-01 9.574954e-01 9.385503e-01 
       8.776899e-01 9.604303e-01 9.397493e-01 
       8.846586e-01 9.632664e-01 9.411828e-01 
       8.917232e-01 9.660086e-01 9.428583e-01 
       8.988707e-01 9.686622e-01 9.447827e-01 
       9.060882e-01 9.712326e-01 9.469623e-01 
       9.133623e-01 9.737257e-01 9.494024e-01 
       9.206794e-01 9.761475e-01 9.521080e-01 
       9.280258e-01 9.785041e-01 9.550830e-01 
       9.353875e-01 9.808020e-01 9.583309e-01 
       9.427505e-01 9.830480e-01 9.618541e-01 
       9.501006e-01 9.852488e-01 9.656544e-01 
       9.574237e-01 9.874116e-01 9.697328e-01 
       9.647055e-01 9.895434e-01 9.740896e-01 
       9.719318e-01 9.916517e-01 9.787242e-01 
       9.790885e-01 9.937439e-01 9.836351e-01 
       9.861613e-01 9.958275e-01 9.888201e-01 
       9.931365e-01 9.979103e-01 9.942764e-01 
       1 1 1 ] ; 