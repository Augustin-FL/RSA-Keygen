#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>


int algo_euclide_etendu(long a_long,long b_long, long *pgcd,long *u_long);
long corrige_u_besout(int uint, int bint);
void menu_crypter_decrypter(long e,long n,long d);

int main(int argc, char **argv)
{
	long p_long,q_long,phi_long,e,pgcd,d,n_long;


	printf("Choisissez deux nombres P et Q premiers entre eux : \np= ");

	scanf("%ld",&p_long);

	printf("q= ");
	scanf("%ld",&q_long);
	
	algo_euclide_etendu(p_long,q_long,&pgcd,&d);

	if(pgcd==1)
	{
		mpz_t p, q,n,phi;
		mpz_init(phi);
		mpz_init(n);
		
		mpz_init_set_si(p,p_long);
		mpz_init_set_si(q,q_long);
		
		printf("calcul de N...\n");
		mpz_mul(n,p,q);// N=P*Q
		
		mpz_sub_ui(p,p,1);
		mpz_sub_ui(q,q,1);
		
		printf("calcul de phi...\n");//phi = (Q-1)(P-1)
		mpz_mul(phi,p,q);
		
		mpz_add_ui(p,p,1);
		mpz_add_ui(q,q,1);
		
		phi_long=mpz_get_si(phi);
		n_long=mpz_get_si(n);
		
		printf("Choisissez un nombre premier avec phi=%ld\ne= ",phi_long);
		scanf("%ld",&e);

        printf("calcul de D...\n");//d est u tel que phi*u+e*v=1 (il y a forcément une infinité de solution car e et phi sont premiers)
		algo_euclide_etendu(phi_long,e,&pgcd,&d);

		if(pgcd == 1)
		{
			printf("Correction des coeficients de bezouts de D....\n");//les coefs de bezout ne peuvent pas être negatifs ni supérieur au modulo.
			d=corrige_u_besout(d,phi_long);
		
			printf("\n\n\nRESULTATS DU CALCUL DES CLES :\n-------------------\nn=%ld\ne=%ld\n\nd (cle secrete) : %ld\n\n",n_long,e,d);
			printf("N et E sont les cles publiques. D est la cle secrete.\nP et Q ont etes detruits de la memoire, et ne devraient etre enregistre nul part.\n\n");
		
			mpz_clear(p);
			mpz_clear(q);
			p_long=0;
			q_long=0;
		
			menu_crypter_decrypter(e,n_long,d);
		}
		else printf("ERREUR : phi et E doivent etre premiers.\nCe n'est pas le cas ici. (PGCD(%ld,%ld) = %ld)",phi_long,e,pgcd);
				
	}
	else  printf("ERREUR : P et Q doivent etre premiers.\nCe n'est pas le cas ici. (PGCD(%ld,%ld)=%ld)",p_long,q_long,pgcd);
	
	return 0;
}


void menu_crypter_decrypter(long e,long n,long d)
{
	int continuer=1,choix=0,merci_j_et_m=1;
	long x_l;
	mpz_t x;

	while(continuer && merci_j_et_m)
	{
		printf("\n\nMENU PRINCIPAL\n---------------------\n1 : Crypter un chiffre\n2 : Decrypter un chiffre\nAutre : Quitter\n\nVotre choix : ");

		merci_j_et_m=scanf("%d",&choix);

		if(choix==1)
		{
			printf("entrez le chiffre a crypter: ");
			scanf("%ld",&x_l);

			mpz_init_set_ui(x,x_l);// pour crypter : x^e % n
			mpz_pow_ui(x,x,e);
			mpz_mod_ui(x,x,n);

			printf("Chiffre crypte : %ld\n\n\n",mpz_get_ui(x));
		}
		else if(choix==2)
		{
			printf("entrez le chiffre a decrypter: ");
			scanf("%ld",&x_l);

			mpz_init_set_ui(x,x_l);// pour décrypter : x^d %n
			mpz_pow_ui(x,x,d);
			mpz_mod_ui(x,x,n);

			printf("Chiffre decrypte : %ld\n\n\n",mpz_get_ui(x));
        }
		else continuer=0;
	}
}


int algo_euclide_etendu(long a_long,long b_long, long *pgcd,long *u_long)
{
	// --------------------- ALGO EUCLIDE ETENDU
	//tiré depuis wikipedia.

	mpz_t a, b, u, v, x, y, c, d, q, r,s,q1,q2;

	mpz_init_set_si(a,a_long);
	mpz_init_set_si(b,b_long);
	mpz_init(q1);
	mpz_init(q2);
	mpz_init(s);
	mpz_init(c);
	mpz_init(d);
	mpz_init(q);
	mpz_init(r);

	mpz_init_set_si(u,1);
	mpz_init_set_si(v,0);

	mpz_init_set_si(x,0);
	mpz_init_set_si(y,1);


	while (mpz_cmp_si(b,0))//while b>0
	{

		mpz_mod(r,a,b);
		mpz_sub(s,a,r);
		mpz_div(q,s,b);
		mpz_set(c,u);
		mpz_set(d,v);
		mpz_set(u,x);
		mpz_set(v,y);

		mpz_mul(q1,q,x);
		mpz_mul(q2,q,y);
		mpz_sub(x,c,q1);
		mpz_sub(y,d,q2);
		mpz_set(a,b);
		mpz_set(b,r);
	}
	
    *pgcd=mpz_get_si(a);
    *u_long=mpz_get_si(v);

    return 0;
}

long corrige_u_besout(int uint, int bint)
{
   /* RSA impose que 2 < u < b...
	 * donc si u < 2 ou u > b, on calcule l'inverse de u :p
	  on a a*(u - k*b) + b*(v + k*a) = 1 -> il suffit de faire varier k
	 */

    mpz_t u,b,k,kb;

    mpz_init_set_si(u,uint);
    mpz_init_set_si(b,bint);
    mpz_init(kb);


	if(uint < 2)
	{
        mpz_init_set_si(k,2);

		mpz_sub(k,k,u);
		mpz_mul_si(b,b,-1);

		mpz_fdiv_q(k, k, b);// k =(2-u)/-b

		mpz_mul_si(b,b,-1);
        mpz_mul(kb,k,b);
		mpz_sub(u,u,kb);
	}

	if(uint > bint)
	{
		mpz_set_si(k,b);


		mpz_sub(k,k,u);
		mpz_mul_si(b,b,-1);

        mpz_fdiv_q(k, k, b);//k=(b-u)/-b

		mpz_mul_si(b,b,-1);
        mpz_mul(kb,k,b);
		mpz_sub(u,u,kb);
	}

	return mpz_get_si(u);
}





