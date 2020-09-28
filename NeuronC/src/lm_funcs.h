
void lmfit(double (*user_func) (double user_x_point, double *coeff, int model_num),
		                int m_dat, double *xydata, int n_p, double *coeff, double *coeffc);

void lmfit2d(double (*user_func2d) (double user_x_point, double user_y_point, double *coeff, int model_num),
		                int m_dat, double *xydata, int n_p, double *coeff, double *coeffc);

