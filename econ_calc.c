/*Example
Input
 1200
Enter time: 2
Enter rate: 5.4
Output

Simple Interest = 129.600006*/
#include<stdio.h>
int main(){
	double p, t, r;
	printf("Enter principle:\n");
	scanf("%lf ",&p);
	printf("Enter time: ");
	scanf("%lf ",&t );
	printf("Enter rate: ");
	scanf("%lf ", &r);


	double answer = (p*t*r)/100;

	printf("Simple Interest = %lf",answer);

	return 0;




}
