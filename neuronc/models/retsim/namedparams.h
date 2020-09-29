#include <string>
#include <vector>
#include <map>

using namespace std;

class namedparams
{

 public:
  namedparams(const char *fname);

  double get(const char *reg, const char *param);
  double get(string reg, string param);
  double get(const char *reg, const char *param, double def);
  double get(string reg, string param, double def);

  bool   has(char *reg, const char *param);
  bool   has(string reg, string param);

 protected:
  map<string, int> nameindex;
  map< string, vector<double> *> paramvals;
};
