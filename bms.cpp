/**
 * Author: Michal Dostal
 * Login: xdosta51
 * e-mail: xdosta51@stud.fit.vutbr.cz
 * file: bms.cpp
 * project: implement RSC
*/


/*
Implementation inspired from source in assigment viz.: https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders
This implementation doesnt work for our assigment, there needs to be made some changes...

*/


//Includes libraries
#include <cmath>
#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <bitset>
#include <iterator>
#include <algorithm>

/*
From assigment used prefilled tables
Preffiled tables: viz. https://github.com/lrq3000/unireedsolomon/blob/master/unireedsolomon/ff.py
*/


// GF exptable prefilled
int gf2int_exptable [256] = {1, 3, 5, 15, 17, 51, 85, 255, 26, 46, 114, 150, 161, 248, 19,
        53, 95, 225, 56, 72, 216, 115, 149, 164, 247, 2, 6, 10, 30, 34,
        102, 170, 229, 52, 92, 228, 55, 89, 235, 38, 106, 190, 217, 112,
        144, 171, 230, 49, 83, 245, 4, 12, 20, 60, 68, 204, 79, 209, 104,
        184, 211, 110, 178, 205, 76, 212, 103, 169, 224, 59, 77, 215, 98,
        166, 241, 8, 24, 40, 120, 136, 131, 158, 185, 208, 107, 189, 220,
        127, 129, 152, 179, 206, 73, 219, 118, 154, 181, 196, 87, 249, 16,
        48, 80, 240, 11, 29, 39, 105, 187, 214, 97, 163, 254, 25, 43, 125,
        135, 146, 173, 236, 47, 113, 147, 174, 233, 32, 96, 160, 251, 22,
        58, 78, 210, 109, 183, 194, 93, 231, 50, 86, 250, 21, 63, 65, 195,
        94, 226, 61, 71, 201, 64, 192, 91, 237, 44, 116, 156, 191, 218,
        117, 159, 186, 213, 100, 172, 239, 42, 126, 130, 157, 188, 223,
        122, 142, 137, 128, 155, 182, 193, 88, 232, 35, 101, 175, 234, 37,
        111, 177, 200, 67, 197, 84, 252, 31, 33, 99, 165, 244, 7, 9, 27,
        45, 119, 153, 176, 203, 70, 202, 69, 207, 74, 222, 121, 139, 134,
        145, 168, 227, 62, 66, 198, 81, 243, 14, 18, 54, 90, 238, 41, 123,
        141, 140, 143, 138, 133, 148, 167, 242, 13, 23, 57, 75, 221, 124,
        132, 151, 162, 253, 28, 36, 108, 180, 199, 82, 246, 1};

// GF Logtable prefilled
int GF2int_logtable [256] = {-1, 0, 25, 1, 50, 2, 26, 198, 75, 199, 27, 104, 51, 238, 223, 
        3, 100, 4, 224, 14, 52, 141, 129, 239, 76, 113, 8, 200, 248, 105,
        28, 193, 125, 194, 29, 181, 249, 185, 39, 106, 77, 228, 166, 114,
        154, 201, 9, 120, 101, 47, 138, 5, 33, 15, 225, 36, 18, 240, 130,
        69, 53, 147, 218, 142, 150, 143, 219, 189, 54, 208, 206, 148, 19,
        92, 210, 241, 64, 70, 131, 56, 102, 221, 253, 48, 191, 6, 139, 98,
        179, 37, 226, 152, 34, 136, 145, 16, 126, 110, 72, 195, 163, 182,
        30, 66, 58, 107, 40, 84, 250, 133, 61, 186, 43, 121, 10, 21, 155,
        159, 94, 202, 78, 212, 172, 229, 243, 115, 167, 87, 175, 88, 168,
        80, 244, 234, 214, 116, 79, 174, 233, 213, 231, 230, 173, 232, 44,
        215, 117, 122, 235, 22, 11, 245, 89, 203, 95, 176, 156, 169, 81,
        160, 127, 12, 246, 111, 23, 196, 73, 236, 216, 67, 31, 45, 164,
        118, 123, 183, 204, 187, 62, 90, 251, 96, 177, 134, 59, 82, 161,
        108, 170, 85, 41, 157, 151, 178, 135, 144, 97, 190, 220, 252, 188,
        149, 207, 205, 55, 63, 91, 209, 83, 57, 132, 60, 65, 162, 109, 71,
        20, 42, 158, 93, 86, 242, 211, 171, 68, 17, 146, 217, 35, 32, 46,
        137, 180, 124, 184, 38, 119, 153, 227, 165, 103, 74, 237, 222, 197,
        49, 254, 24, 13, 99, 140, 128, 192, 247, 112, 7};



// Function Galois fields sub
// parameters are int x, y
// it outputs int 
int gf_sub(int x, int y) {
    return x ^ y;
}

// GF exp multiplication, check oversize
// parameters are int x, y
// it outputs int 
int gf_mul(int x, int y) {
    
    if (x == 0 || y == 0) 
        return 0;
    
    return gf2int_exptable[(GF2int_logtable[x] + GF2int_logtable[y]) % 255];
}

// GF add function
// parameters are int x, y
// it outputs int 
int gf_add(int x, int y) {
    return pow(x,y);
}

// Multiply polynomials
// Parametres are 2 polynoms which are stored in vector of int
// It returns vector of ints
std::vector<int> gf_poly_mul(std::vector<int> p, std::vector<int> q) {
    // fill array with zeros
    std::vector<int> temp_arr ((p.size() + q.size()) -1);
    for (long unsigned int i = 0; i < ((p.size() + q.size()) -1); i++) {
        temp_arr[i] = 0;
    }
    // calculate temp_arr
    for (long unsigned int j = 0; j < (q.size()); j++) {
       for (long unsigned int i = 0; i < (p.size()); i++) {
            
            temp_arr[i + j] ^= gf_mul(p[i], q[j]);
        } 
    }

    return temp_arr;
}

// Galois fields div function
// parameters are int x, y
// it returns y
int gf_div(int x, int y) {
    if (y == 0) {
        std::cerr << "Zero division" << "\n";
        exit(0);
    }
    if (x == 0) {
        return (0);
    }
    return gf2int_exptable[(GF2int_logtable[x] + 255 - GF2int_logtable[y]) % 255]; // modulo needs to be added
}

//GF power function
// arguments are int x and int power
// it returns int 
int gf_pow(int x, int power) {
    // negative indexing repair
    if ((GF2int_logtable[x] * power) < 0)
        return gf2int_exptable[255 + ((GF2int_logtable[x] * power) % 255)];
    return gf2int_exptable[(GF2int_logtable[x] * power) % 255];
}

// galois field invers
// parameters are int x
// it returns int
int gf_inverse(int x) {
    return gf2int_exptable[255 - GF2int_logtable[x]];
}


//Galois field polynomial scale
//parameters are polynom p, whic is vector of int, and int x
//it returns polynom
std::vector<int> gf_poly_scale(std::vector<int> p, int x) {
    std::vector<int> temp_arr(p.size());
    //fill array with zeros.
    for (long unsigned int i = 0; i < p.size(); i++) {
        temp_arr[i] = 0;
    }
    for (long unsigned int i = 0; i < p.size(); i++) {
        temp_arr[i] = gf_mul(p[i], x);
    }
    return temp_arr;


}


//Galois field polynomial add
//parameters are 2 polynoms
//output is one polynom which is result off add of 2 polynoms from parameters
std::vector<int> gf_poly_add(std::vector<int> p, std::vector<int> q) {
    int max = 0;
    if (p.size() > q.size())
        max = p.size();
    else 
        max = q.size();
    // fill array with zeros
    std::vector<int> temp_arr(max);
    for (int i = 0; i < max; i++) {
        temp_arr[i] = 0;
    }
    // calculate array
    for (long unsigned int i =0; i < p.size(); i++) {
        temp_arr[i+max - p.size()] = p[i];
    }
    for (long unsigned int i =0; i < q.size(); i++) {
        temp_arr[i+max - q.size()] ^= q[i];
    }
    
    return temp_arr;

}

// Polynom generator for value gen=3 and fcr=1
// This function needs to have change value in gf_pow for generator
// it takes one argument which is size of generator
// it returns polynom
std::vector<int> rs_generator_poly(int nsym) {
    // calculate generator
    std::vector<int> g (nsym);

    g = {1};

    for (int i = 0; i < nsym; i++) {
        g = gf_poly_mul(g, {1, gf_pow(3, i+1)});
    }

    return g;
}

// GF Polynomial eval function
// arguments are polynom and int x
// int returns int
int gf_poly_eval(std::vector<int> poly, int x) {
    
    int y = poly[0];
    
    for (long unsigned int i = 1; i < poly.size(); i++) {
        y = gf_mul(y,x) ^ poly[i];
    }
    
    return y;
}

// Polynomial div function
// This function is for encoding only because it works with string
// input is string, generator and nsym whic is computed size of n and k
// it returns remainder, which is the string added to source string.
std::string gf_poly_div (std::string dividend, std::vector<int> generator, int nsym) {
    
    std::vector<int> msg_out (dividend.size() + nsym);

    for (long unsigned int i = 0; i < dividend.size(); i++) {
        msg_out[i] = int(dividend[i]);
    }

    for (long unsigned int i = 0; i < dividend.size() - (generator.size() - 1); i++) {
        int coef = msg_out[i];
        
        if (coef != 0) {
            for (long unsigned int j = 1; j < generator.size(); j++) {
                if (generator[j] != 0) {
                    
                    
                     msg_out[i + j] ^= gf_mul(generator[j], coef);
                    
                }
            }
        }

    }
    

    std::string msg_out_final;

    for (long unsigned int i = 0; i < dividend.size(); i++) {
        msg_out_final += msg_out[i];
        
    }
    // return only remainder
    return msg_out_final.substr(msg_out_final.length() - (nsym));
}


//rs find errata locator
// This function find err locator
// parameters is position where errors are located
// it returns vector of ints.
// gf pow must be changed to 3 because of gen=3
std::vector<int> rs_find_errata_locator(std::vector<int> e_pos) {
    std::vector<int> e_loc = {1};

    for (long unsigned int i = 0; i < e_pos.size(); i++) {
        std::vector<int> e_loc_temp = gf_poly_mul(e_loc, gf_poly_add(std::vector<int> {1}, std::vector<int> {gf_pow(3, e_pos[i]), 0}));

        
        e_loc.resize(e_loc_temp.size());
        e_loc = e_loc_temp;;
    }
    

    return e_loc;
}



// Galois field polynomial div decode
// This function is only for decoding
// It works the same as gf_poly_div but with vector of ints not with string
// and it doesnt need param nsym
// it returns remainder
std::vector<int> gf_poly_div_decode(std::vector<int> dividend, std::vector<int> divisor) {
    std::vector<int> msg_out = dividend;
    for (long unsigned int i = 0; i < dividend.size() - (divisor.size() - 1); i++) {
        int coef = msg_out[i];
        if (coef != 0) {
            for (long unsigned int j = 1; j < divisor.size(); j++) {
                if ((divisor[j]) != 0)
                    msg_out[i+j] ^= gf_mul(divisor[j], coef);
            }
        }
    }
    int separator = -(divisor.size() - 1);

    std::vector<int> y(msg_out.end() + separator, msg_out.end());

    return y;

}

// Reed solomon find error evaluator
// inputs are syndroms, err_locations, and nsym
// it returns vector of ints.
std::vector<int> rs_find_error_evaluator(std::vector<int> synd, std::vector<int> err_loc, int nsym) {
    
    std::vector<int> first_arg = gf_poly_mul(synd, err_loc);
    std::vector<int> second_arg (nsym+2);

    second_arg[0] = 1;
    // fill array with zeros at start is 0
    for (int i = 0; i < nsym+1; i++)
        second_arg[i+1] = 0;

    std::vector<int> remainder = gf_poly_div_decode(first_arg, second_arg);
    return remainder;

}



//Function to correct errors in message
// This function outputs corrected message
// It has arguments msg_in which is message with errors, list of syndromes, and err_positions.
// gf pow must be change to 3 because of gen =3
// y2 gf pow second argument changed to 0 because of assigment
std::vector<int> rs_correct_errata(std::vector<int> msg_in, std::vector<int> synd, std::vector<int> err_pos) {
    std::vector<int> coef_pos (err_pos.size());

    
    
    // inverse coef_pos
    for (long unsigned int i = 0; i < err_pos.size(); i++) {
        
        coef_pos[coef_pos.size() - i - 1] = msg_in.size() - 1 - i;
    }

    
    //find err_locat0r
    std::vector<int> err_loc = rs_find_errata_locator(coef_pos);
    //reverse syndroms list
    std::reverse(synd.begin(), synd.end());
    // find_errors and revers the vector
    std::vector<int> err_eval = rs_find_error_evaluator(synd, err_loc, err_loc.size()-1);
    std::reverse(err_eval.begin(), err_eval.end());

    
   


    std::vector<int> X (coef_pos.size());
    // calculate X
    for (long unsigned int i = 0; i < coef_pos.size(); i++) {
        int l = 255 - coef_pos[i];
        
        
        X[i] =  gf_pow(3, -l);
        
        
    }
   
    // fill E with zeros
    std::vector<int> E (msg_in.size());
    for (long unsigned int i =0; i < msg_in.size(); i++) {
        E[i] = 0;
    }

    int Xlength = X.size();

    for (long unsigned int i = 0; i < X.size(); i++) {
        
        int Xi = X[i];
        int Xi_inv = gf_inverse(Xi);
        
        
        std::vector<int> err_loc_prime_tmp;

        for (int j = 0; j < Xlength; j++) {
            if (j != int(i)) {
                err_loc_prime_tmp.push_back(gf_sub(1, gf_mul(Xi_inv, X[j])));
                
            }
        }

        int err_loc_prime = 1;
        for (long unsigned int iter = 0; iter < err_loc_prime_tmp.size(); iter++) {
            err_loc_prime = gf_mul(err_loc_prime, err_loc_prime_tmp[iter]);
            
            
        }
        // reverse list
        std::reverse(err_eval.begin(), err_eval.end());
        
        int y = gf_poly_eval(err_eval, Xi_inv);
       
        int y2 = gf_mul(gf_pow(Xi, 0), y);
        

        if (err_loc_prime == 0) {
            std::cerr << "err";
            exit(0);
        }
        // reverse list to back to normal
        std::reverse(err_eval.begin(), err_eval.end());
        // calculate repaired symbol
        int magnitude = gf_div(y2, err_loc_prime);
        E[err_pos[i]] = magnitude;
    }

    msg_in = gf_poly_add(msg_in, E);
    return(msg_in);

}

//Functions to find locations with errors
// Inputs are syndromes and nsym
// It outputs vector of ints with err_location
std::vector <int> rs_find_error_locator(std::vector<int> synd, int nsym) {
    std::vector<int> err_loc = {1};
    std::vector<int> old_loc = {1};

    long unsigned int synd_shift = synd.size() - nsym;

    for (int i = 0; i < nsym; i++) {
        int K = i + synd_shift;
        int delta = synd[K];
        for (int j = 1; j < int(err_loc.size()); j++) {
            delta ^= gf_mul(err_loc[err_loc.size()-(j+1)], synd[(K - j)]);
        }

        old_loc.resize(old_loc.size() + 1);
        old_loc[old_loc.size()] = 0;

        if (delta != 0) {
            if (old_loc.size() > err_loc.size()) {
                std::vector<int> temp_new_loc = gf_poly_scale(old_loc, delta);
                
                std::vector<int> temp_old_loc = gf_poly_scale(err_loc, gf_inverse(delta));
                old_loc.resize(temp_old_loc.size());
                old_loc = temp_old_loc;
                
                

                err_loc.resize(temp_new_loc.size());
                err_loc = temp_new_loc;

            }
            std::vector<int> temp_err_loc = gf_poly_add(err_loc, gf_poly_scale(old_loc, delta));
            err_loc.resize(temp_err_loc.size());
            err_loc = temp_err_loc;
            
            
        } 

        for (long unsigned int i =0; i < err_loc.size(); i++) {
        //std::cout << "ERRLOC" << err_loc[i] << " \n";
    }
    }

    while (err_loc.size() && (err_loc[0] == 0)) 
        err_loc.erase(err_loc.begin());
    long unsigned int errs = (err_loc.size() - 1);
    if (int((errs * 2)) > nsym) {
        std::cerr << "Too many errs\n";
        exit(0);
    }

    return err_loc;
}

// Function to find errors
// This functions find errors in source string
// it does by the line if (gf_poly_eval(err_loc, gf_pow(3, i)) == 0) this line has to be edited because of gen=3
// if this condition is zero it has error 
// It returns positions with errors.
std::vector<int> rs_find_errors(std::vector<int> err_loc, int nmess) {
    int errs = err_loc.size() - 1;
    std::vector<int> err_pos;

    for (int i = 0; i < nmess; i++) {
        if (gf_poly_eval(err_loc, gf_pow(3, i)) == 0) { // if this condition is true there is error
            err_pos.resize(err_pos.size()+1);
            err_pos[err_pos.size()-1] = nmess - 1 - i;
            
        }
    }
    if (int(err_pos.size()) != errs) { 
        std::cerr << "Too many errs\n";
        exit(0);
    }
    return err_pos;

}

//Encoding function
// This functions take source string to encode and nsym.
// It prints the final message on output in binary 
// and returns 0
int encode(std::string string_to_encode, int nsym) {
    
    std::vector<int> gen (nsym);
    gen = rs_generator_poly(nsym);
    // fill array with zeros
    std::string temp_2;
    for (int i = 0; i < (nsym); i++) {
        temp_2 += char(0);
    }
    // calculate encoded remainder
    std::string remainder = gf_poly_div(string_to_encode + temp_2, gen, nsym);
    // make final string
    std::string msg_out = string_to_encode + remainder;
    // print it in binary
    for (long unsigned int i = 0; i < msg_out.size(); i++) {
        std::bitset<8> x(msg_out[i]);
        std::cout << x;
    }
    // print newline at end
    std::cout << "\n";
    
    return 0;
}
//Calculate syndromes..
// This function returns list with zeros if no error is in message
// or returns non zero list if message has errors.
// Inputs are msg which is source message and nsym, which is computated from n and k
// gf pow is 3, i+1 becuase of gen=3 and fcr=1
std::vector<int> rs_calc_syndromes(std::vector<int> msg, int nsym) {
    
    std::vector<int> temp_arr (nsym);

    for (int i = 0; i < nsym; i++) {
        temp_arr[i] = 0;
    }

    for (int i =0; i < nsym; i++) {
        

        temp_arr[i] = gf_poly_eval(msg, gf_pow(3,i+1)); // gen=3 , fcr=1
    }
    // first syndrom is zero
    // 
    std::vector<int> temp_arr_out (nsym+1);
    for (int i = 0; i < nsym+1; i++) {
        if (i == 0) 
            temp_arr_out[i] = 0;
        else 
            temp_arr_out[i] = temp_arr[i-1];
    }
    return temp_arr_out;

}


// Function for decoding..
// this functions decode message encoded in reed solomon code with generator=3
// it has args with nsym, source encoded string and size of final string.
// it returns 0 if everthings its ok and print decoded message, else return 1 and print err
int decode(int decode_nsym_arg, std::string string_to_decode, int size_of_decoded_string) {
    // declare variable for string to decode from binary
    std::vector<int> array_of_msg_to_decode (string_to_decode.size()/8);
    int save_to_array = 0;

    // this for makes from binary representation list with ints from 1 byte
    for (long unsigned int i = 0; i < string_to_decode.size(); i=i) {
        std::string temp_string = string_to_decode.substr(i,8);
        std::bitset<8> b(temp_string); 
        unsigned char c = ( b.to_ulong() & 0xFF);
       
        array_of_msg_to_decode[save_to_array] = int(c);
        save_to_array++;
        i = i+8;
    }
    // calculate syndromes
    std::vector<int> syndromes = rs_calc_syndromes(array_of_msg_to_decode, decode_nsym_arg);

    int string_ok = 0;

    for (long unsigned int i = 0; i < syndromes.size(); i++) {
        if (syndromes[i] != 0)
            string_ok = 1;
    }
    // if max syndrome was 0 the message is ok we cant print the msg without remainder and exit program
    if (!string_ok) {
        for (int i = 0; i < size_of_decoded_string; i++) 
            std::cout << char(array_of_msg_to_decode[i]);
        std::cout << "\n";
        return(0);
    }

    
    // remove 0 from vector of syndromes
    std::vector <int> new_syndromes(syndromes.size()-1);
    for (long unsigned int i = 0; i < new_syndromes.size(); i++) {
        new_syndromes[i] = syndromes[i+1];
        
    }
    // find errors
    std::vector<int> find_errors = rs_find_error_locator(new_syndromes, decode_nsym_arg);
    // reverse vector with errors
    std::reverse(find_errors.begin(), find_errors.end());
    // find err_positions
    std::vector<int> err_pos = rs_find_errors(find_errors, array_of_msg_to_decode.size());
    // correct the message.
    std::vector<int> corrected = rs_correct_errata(array_of_msg_to_decode, syndromes, err_pos);
    // check the correct message for syndromes
    std::vector<int> calc_new_syndromes = rs_calc_syndromes(corrected, decode_nsym_arg);

    int without_errs = 0;

    for (long unsigned int i = 0; i < calc_new_syndromes.size(); i++) {
        if (calc_new_syndromes[i] != 0) {
            without_errs = 1;
        }
    }

    // if corrected message had max syndrome 0 its ok and we can print it and exit program.
    if (!without_errs) {
        for (int i = 0; i < size_of_decoded_string; i++) 
            std::cout << char(corrected[i]);
    }
    else {
        std::cerr << "It contains errors \n";
        return(1);
    }
    std::cout << "\n";

    return 0;

}

// main function, parse arguments
int main(int argc, char* argv[]) {

    if (argc > 1 && argc < 9) {
        
        if (argv[1] == std::string("-e")) {
            
            std::string string_to_encode;
            int size_of_encoded_string;
            int size_of_string_to_encode;

            if (argv[2] == std::string("-n")) {

                size_of_encoded_string = atoi(argv[3]);
                string_to_encode = argv[5];

            }

            else {

                size_of_encoded_string = atoi(argv[5]);
                string_to_encode = argv[3];

            }
            size_of_string_to_encode = string_to_encode.size();
            int nsym = size_of_encoded_string - size_of_string_to_encode;

            encode(string_to_encode, nsym);
            
            return(0);
        }

        else if (argv[1] == std::string("-d")) {
            
            std::string string_to_decode;
            int size_of_decoded_string;
            int size_of_string_to_decode;

            if (argc != 8) {
                std::cerr << "Spatny pocet argumentu\n";
                return(1);
            }
            
            for (int i = 2; i < 8; i++) {
                if (argv[i] == std::string("-n")) {
                    size_of_string_to_decode = atoi(argv[i+1]);
                    i++;
                }
                else if (argv[i] == std::string("-k")) {
                    size_of_decoded_string = atoi(argv[i+1]);
                    i++;
                }
                else if (argv[i] == std::string("-m")) {
                    string_to_decode = argv[i+1];
                    i++;
                }
            }

            int decode_nsym = size_of_string_to_decode - size_of_decoded_string;

            decode(decode_nsym, string_to_decode, size_of_decoded_string);
            
            return(0);
        }
    }

    else if (argc >= 9) {
        std::cerr << "Moc velky pocet argumentu\n";
        return(1);
    }
    
    return 0;
}