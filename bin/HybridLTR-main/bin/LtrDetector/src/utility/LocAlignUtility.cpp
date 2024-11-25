#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>

#include "../utility/LocAlign.h"
using namespace std;
using namespace utility;


int main(int argc, char *argv[])
{

    ifstream ifs;

    ifs.open(argv[1], ifstream::in);
    cout << "FILE OPENED" << endl;
    char c = ifs.get();

    if (c == '>')
    {

        while (c != '\n')
        {
            c = ifs.get();
        }
    }

    string string1 = "";

    while (ifs.good())
    {

        if (c != '\n')
        {
            string1 += c;
        }
        c = ifs.get();
    }

    ifs.close();

    ifstream ifs2;

    ifs2.open(argv[2], ifstream::in);

    c = ifs2.get();

    if (c == '>')
    {

        while (c != '\n')
        {
            c = ifs2.get();
        }
    }

    string string2 = "";

    while (ifs2.good())
    {


        if (c != '\n')
        {
            string2 += c;
        }
        c = ifs2.get();
    }

    ifs2.close();

    std::transform(string1.begin(), string1.end(), string1.begin(), ::toupper);
    std::transform(string2.begin(), string2.end(), string2.begin(), ::toupper);

    //cout <<string1<<endl;
    //cout <<string2<<endl;

    LocAlign *align = new LocAlign(string1.c_str(), 0, string1.size() - 1, string2.c_str(), 0, string2.size() - 1, 2, -3, 5, 2,false);
    cout << "FINAL SCORE =" << align->getScore() << endl;
    align->printAlignment();
    cout << "ALIGNMENT LENGTH =" << align->getLength() << endl;
    cout << "IDENTITY" << align->getIdentity() << endl;
}
