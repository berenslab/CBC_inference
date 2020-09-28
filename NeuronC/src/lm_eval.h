typedef struct {
    double *user_t;
    double *user_y;
    union {
        double (*user_func)   (double user_t_point, double *par, int model_num);
        double (*user_func2d) (double user_t_pointx, double user_t_pointy, double *par, int model_num);
        struct Symbol *user_funcp;
        struct Symbol *user_funcp2d;
    };
    int user_n_p;
} lm_data_type;

void lm_evaluate_default(double *par, int m_dat, double *fvec,
			 void *data, int *info, int model_num);

void lm_evaluate_2d(double *par, int m_dat, double *fvec,
			 void *data, int *info, int model_num);

void lm_evaluate_interp(double *par, int m_dat, double *fvec,
			 void *data, int *info, int model_num);

void lm_evaluate_2d_interp(double *par, int m_dat, double *fvec,
			 void *data, int *info, int model_num);

void lm_print_default(int n_par, double *par, int m_dat, double *fvec,
		      void *data, int iflag, int iter, int nfev);

void lm_print_2d(int n_par, double *par, int m_dat, double *fvec,
		      void *data, int iflag, int iter, int nfev);

void lm_print_interp(int n_par, double *par, int m_dat, double *fvec,
		      void *data, int iflag, int iter, int nfev);

void lm_print_2d_interp(int n_par, double *par, int m_dat, double *fvec,
		        void *data, int iflag, int iter, int nfev);

