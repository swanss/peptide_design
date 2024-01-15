#ifndef mstsystem_exts_h
#define mstsystem_exts_h

#include "mstsystem.h"

using namespace std;
using namespace MST;

class MstSystemExtension: public MstSys {
public:
    // Returns the file name without directories (including extension)
    static string fileName(const string& path) {
        if (path.find_last_of("/") == string::npos) return path;
        else return path.substr(path.find_last_of("/") + 1);
    }
    
    static string join(const string& path1, const string& path2) {
        if (stringEndsWith(path1, "/") && stringStartsWith(path2, "/")) {
            return path1 + path2.substr(1, path2.length() - 1);
        } else if (stringEndsWith(path1, "/") || stringStartsWith(path2, "/")) {
            return path1 + path2;
        }
        return path1 + "/" + path2;
    }

    /**
     Returns path1 but relative to basePath. Requires that path1 be a child of
     basePath in the directory tree. For example, /a/b/c/ is a child of /a/, but
     /a/b/d is not a child of /a/b/c.
     
     @param path1 the absolute path
     @param basePath the path to be relative to
     @return path1 but relative to basePath
     */
    static string relativePath(const string& path1, const string& basePath) {
        string base = basePath;
        
        if (!stringEndsWith(base, "/"))
            base = base + "/";
        MstUtils::assertCond(stringStartsWith(path1, base), "path must begin with base path");
        
        return path1.substr(base.length(), path1.length() - base.length());
    }

private:
    static bool stringEndsWith(std::string const &fullString, std::string const &ending) {
        if (fullString.length() >= ending.length()) {
            return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
        } else {
            return false;
        }
    }

    static bool stringStartsWith(std::string const &fullString, std::string const &starting) {
        if (fullString.length() >= starting.length()) {
            return (0 == fullString.compare (0, starting.length(), starting));
        } else {
            return false;
        }
    }
};

#endif /* mstsystem_exts_h */
