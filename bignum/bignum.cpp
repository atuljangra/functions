#include<iostream>
#include<cstdio>
#include<cstdlib>

// #define NDEBUG
#include<cassert>

#define MAXLEN 10
class bignum{
    char data[MAXLEN];
    int length;
  public:
  
  void make_bignum (const char* s, bignum *bn) {
    for (bn->length = 0; s[bn->length]&&(bn->length < MAXLEN); bn->length++)
      bn->data[bn->length] = s[bn->length];
  }
  
  void make_bignum (long n, bignum *bn) {
    for (bn->length = 0; n && (bn->length < MAXLEN); bn->length++, n/=10)
      bn->data[bn->length] = n%10 + '0'; 
    char temp;
    for (int i = 0, j = bn->length - 1; i < j; i++ , j--) {
      temp = bn->data[j];
      bn->data[j] = bn->data[i];
      bn->data[i] = temp;
      }
  }
  
  void copy (bignum source, bignum *dest) {
    dest->length = source.length;
    for (int i = 0; i < source.length; i++)
      dest->data[i] = source.data[i]; 
  }
  
  void print_num(bignum bn) {
    for (int i = 0; i < bn.length; i++)
      std::cout <<bn.data[i];
    std::cout << std::endl;
  }
  
  void add (bignum one, bignum two, bignum *result) {
    int i = one.length - 1, j = two.length -1, sum = 0, carry = 0, last = MAXLEN - 1;
    char d[MAXLEN];
    for (; i >= 0 && j >= 0;i--, j--, last--) {
      sum = int(one.data[i]) + int(two.data[j]) + carry - 2*int('0') ;
      carry = sum/10;
      d[last] = sum%10 +'0';
    }
    while (i >= 0) {
      sum = int(one.data[i--]) - int('0') + carry;
      carry = sum/10;
      d[last--] = sum%10 + '0';
    }
    while (j >= 0) {
      sum = int(two.data[j--]) - int('0') + carry;
      carry = sum/10;
      d[last--] = sum%10 + '0';
    }
    while(carry > 0) {
      d[last--] = carry%10  + '0';
      carry/=10;
    }
    result->length = MAXLEN - 1 - last;
    for (int i = 0; i< result->length; i++) {
      result->data[i] = d[i+last+1];
    }
  }

  bool isGreater (bignum big, bignum small) {//true if big >= small
    int i = big.length, j = small.length;
    if (i!= j)
      return (i > j);
    else {
      assert(i == j);
      int k = 0;
      while( k < i) {
        if (big.data[k] != small.data[k])
          return (big.data[k] > small.data[k]);
        assert(big.data[k] == small.data[k]);
        k++;
      }
    }
    //return true if both are equal
    return true;
  }

  bool sub (bignum one, bignum two, bignum* result) {
    //Return true if one is greater than two otherwise false
    bool temp_bool = isGreater(one, two);
    if (!temp_bool) {
      sub(two, one, result);
      return false;
    }
    assert(isGreater(one, two));
    int diff = 0, borrow = 0, i = one.length - 1, j = two.length - 1, last = MAXLEN - 1;
    char temp[MAXLEN];

    while (j >= 0) {
      diff = int(one.data[i--]) - int(two.data[j--]) + borrow;
      borrow = 0;
      while (diff < 0) {
        borrow--;
        diff += 10;
      }
      temp[last--] = diff + '0';
    }
    while (i >=0) {
      diff = int(one.data[i--]) - int('0') + borrow;
      borrow = 0;
      while(diff < 0) {
        borrow--;
        diff += 10;
        }
      temp[last--] = diff + '0';
    }
    /*working fine, need to check case of borrow>0??
    * removing all the preceding zeros except the case 
    * when answer is 0, thus MAXLEN - 2.
    */ 
    while ((last < MAXLEN-2) && (temp[last+1]) == '0') {
      last++;
    }
    
    result->length = MAXLEN - 1 - last;
    for (int i = 0; i< result->length; i++) {
      result->data[i] = temp[i+last+1];
    }
     return true;
  }

  void multiply (bignum one, bignum two, bignum* result) {
    //Maximum possible length of result
    result->length = one.length + two.length;  
    if (result->length > MAXLEN) {
      std::cout << "Error: Overflow while multiplication" << std::endl;
      return; 
    }
  }
};

bool isValid(const char s[]) {
    int i = 0;
    while (s[i]) {
      //Comparing ascii values
      if(s[i] < 48 || s[i] > 57)
        return false;
      i++;
    }
    return true;
  }

int main(){
  bignum one,two,sum, diff;
  one.make_bignum("99999", &one);
  two.make_bignum("99999", &two);
  sum.add(one, two, &sum);
  diff.sub(one, two, &diff);
  one.print_num(one);
  two.print_num(two);
  sum.print_num(sum);
  diff.print_num(diff);
}