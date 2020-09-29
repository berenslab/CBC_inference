
#include <string.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

#include "namedparams.h"
#include "nc.h"

double getvarval(const char *name); //defined in ncsetvar.cc


/* 
   Break a string into a list of string tokens separated
   by the specified delimiters
 */
vector<string> *tokenize(string str, const string &delims = " ")
{
	string::size_type lastPos = str.find_first_not_of(delims, 0);
	string::size_type pos     = str.find_first_of(delims, lastPos);

	vector<string> *tokens = new vector<string>();

	while ((pos != string::npos) || (lastPos != string::npos)) {
		tokens->push_back(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delims, pos);
		pos = str.find_first_of(delims, lastPos);
	} 

	return tokens;
}

/*
  Makes a string all upper case
 */
string string2upper(string str)
{
	char *str2 = new char[str.size()+1];
	for (int k = 0; k < str.size(); k++) str2[k] = toupper((int) str.at(k));
	str2[str.size()] = '\0';
	return string(str2);
}

/*
  Makes a string all lower case
 */
string string2lower(string str)
{
	char *str2 = new char[str.size()+1];
	for (int k = 0; k < str.size(); k++) str2[k] = tolower((int) str.at(k));
	str2[str.size()] = '\0';
	return string(str2);
}


namedparams::namedparams(const char *fname)
{
	ifstream in(fname);
	char *buff = new char[8192];

	int k;
	vector<string> *toks, *pvals;  
	string pname, val;
	double dval;
	char *vstr;
	Symbol *sym;
	while (!in.eof()) {
		in.getline(buff, 8192);
		if (strlen(buff) > 0) {
			if (buff[0] != '#') {
				toks = tokenize(buff, ":");
				if (toks->size() > 1) {
					pname = string((*toks)[0]);
					pvals = tokenize((*toks)[1], " \t");

					if (pname == "columns" || pname == "0") {
						for (k = 0; k < pvals->size(); k++) {
							val = string( string2upper( (*pvals)[k] ) );	      
							//fprintf(stderr, "[namedparams] Adding column: '%s', index=%d\n", val.c_str(), k);
							nameindex[val] = k;
						}
					} else {
						vector<double> *cpvals = new vector<double>;
						//needs to be some sort of error checking here...
						for (k = 0; k < pvals->size(); k++) {	      
							val = (*pvals)[k];

							//check whether "val" is a string or number
							if ((val[0] != '.') && (isalpha(val[0]))) {
								//get value of variable corresponding to "val", check command-line first
								vstr = (char *) val.c_str();
								dval = getvarval(vstr);
								if (dval == LARGENUM) {
									//check symbol table second
									sym = lookup(vstr);
									if (sym != NULL) {
										dval = sym->val;		    
									} else {
										//fprintf(stderr, "# Could not find variable '%s', defaulting to zero!\n", vstr);
										dval = 0;
									}
								}

							} else {
								//convert to double, store
								dval = strtod(val.c_str(), NULL);
							}
							cpvals->push_back(dval);
							//fprintf(stderr, "[namedparams] pname='%s', cpvals[%d]=%g\n", pname.c_str(), k, dval);
						} 
						paramvals[pname] = cpvals;   	    
					}
				}
			}
		}
	}
	delete buff;
}

double namedparams::get(string reg, string param)
{
	map<string,int>::iterator iter = nameindex.find(string2upper(reg));
	map<string, vector<double> *>::iterator iter2 = paramvals.find(param);

	if ((iter != nameindex.end()) && (iter2 != paramvals.end())) {
		vector<double> *pvals = paramvals[param];
		if (pvals != NULL) {
			int indx = nameindex[string2upper(reg)];
			//fprintf(stderr, "[regparams.get] reg=%s, param=%s, indx=%d, dval=%g\n",
			//      reg.c_str(), param.c_str(), indx, dval);
			return (*pvals)[indx];    
		}
	}
	return LARGENUM;
}

double namedparams::get(string colname, string param, double def)
{
	double ret = get(colname, param);
	if (ret == LARGENUM) return def;
	return ret;
}

double namedparams::get(const char *colname, const char *param)
{
	return get(string(colname), string(param));
}

double namedparams::get(const char *colname, const char *param, double def)
{
	return get(string(colname), string(param), def);
}

bool namedparams::has(string colname, string param)
{
	map<string,int>::iterator iter = nameindex.find(string2upper(colname));
	map<string, vector<double> *>::iterator iter2 = paramvals.find(param);

	if ((iter != nameindex.end()) && (iter2 != paramvals.end())) return true;
	return false;
}

bool namedparams::has(char *reg, const char *param)
{
	return has(string(reg), string(param));
}
