/*
Copyright (c) 2003-2010, Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include "_kiss_fft_guts.h"
/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */

static void kf_bfly2(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx * Fout2;
    kiss_fft_cpx * tw1 = st->twiddles;
    __declspec(align(16)) kiss_fft_cpx t;
    Fout2 = Fout + m;
    do{
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        C_MUL (t,  *Fout2 , *tw1);
        tw1 += fstride;
        //C_SUB( *Fout2 ,  *Fout , t );
        //C_ADDTO( *Fout ,  t );
        C_ADDSUB(*Fout , *Fout2 , *Fout , t);
        ++Fout2;
        ++Fout;
    }while (--m);
}

static void kf_bfly4(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        const size_t m
        )
{
    kiss_fft_cpx *tw1,*tw2,*tw3;
    __declspec(align(16)) kiss_fft_cpx scratch[6];
    size_t k=m;
    const size_t m2=2*m;
    const size_t m3=3*m;


    tw3 = tw2 = tw1 = st->twiddles;
#if 1
    if(st->inverse) {
        do {
            C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);

            C_MUL(scratch[0],Fout[m] , *tw1 );
            C_MUL(scratch[1],Fout[m2] , *tw2 );
            C_MUL(scratch[2],Fout[m3] , *tw3 );

            //C_SUB( scratch[5] , *Fout, scratch[1] );
            //C_ADDTO(*Fout, scratch[1]);
            C_ADDSUB(*Fout, scratch[5], *Fout, scratch[1]);

            //C_ADD( scratch[3] , scratch[0] , scratch[2] );
            //C_SUB( scratch[4] , scratch[0] , scratch[2] );
            C_ADDSUB(scratch[3] ,scratch[4] ,scratch[0] ,scratch[2]);

            //C_SUB( Fout[m2], *Fout, scratch[3] );
            //C_ADDTO( *Fout , scratch[3] );
            C_ADDSUB(*Fout, Fout[m2], *Fout, scratch[3]);
            tw1 += fstride;
            tw2 += fstride*2;
            tw3 += fstride*3;
            //C_ADDTO( *Fout, scratch[3]);
# if 1
            {__m128d src1, src2, src3, src4, src5, src6;
            __m128d minus1, minus2;
            src1 = _mm_loadr_pd((double*)&(scratch[4].r));
            //src1 = _mm_load_pd((double*)&(scratch[4].r));
            //src1 = _mm_shuffle_pd(src1, src1, _MM_SHUFFLE2(0, 1));
            src2 = _mm_load_pd((double*)&(scratch[5].r));
#  if 1
            //minus1 = _mm_set_pd(0.0, -0.0);
            minus2 = _mm_set_pd(-0.0, 0.0);
            //src3 = _mm_xor_pd(src1, minus1);  // -4.i,  4.r
            src4 = _mm_xor_pd(src1, minus2);  //  4.i, -4.r
#  else
            //minus1 = _mm_set_pd(1.0, -1.0);
            minus2 = _mm_set_pd(-1.0, 1.0);
            //src3 = _mm_mul_pd(src1, minus1);
            src4 = _mm_mul_pd(src1, minus2);
#  endif
            //src5 = _mm_add_pd(src2, src3);
            src5 = _mm_addsub_pd(src2, src1);  // a0-b0, a1+b1
            src6 = _mm_add_pd(src2, src4);
            _mm_store_pd((double*)&(Fout[m].r), src5);
            _mm_store_pd((double*)&(Fout[m3].r), src6);}
# else
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
# endif
            ++Fout;
        }while(--k);
    }else {
        do {
            C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);

            C_MUL(scratch[0],Fout[m] , *tw1 );
            C_MUL(scratch[1],Fout[m2] , *tw2 );
            C_MUL(scratch[2],Fout[m3] , *tw3 );

            //C_SUB( scratch[5] , *Fout, scratch[1] );
            //C_ADDTO(*Fout, scratch[1]);
            C_ADDSUB(*Fout, scratch[5] , *Fout , scratch[1]);

            //C_ADD( scratch[3] , scratch[0] , scratch[2] );
            //C_SUB( scratch[4] , scratch[0] , scratch[2] );
            C_ADDSUB(scratch[3] ,scratch[4] ,scratch[0] ,scratch[2]);

            //C_SUB( Fout[m2], *Fout, scratch[3] );
            //C_ADDTO( *Fout , scratch[3] );
            C_ADDSUB(*Fout, Fout[m2] , *Fout , scratch[3]);
            tw1 += fstride;
            tw2 += fstride*2;
            tw3 += fstride*3;
# if 1
            {__m128d src1, src2, src3, src4, src5, src6;
            __m128d minus1, minus2;
            src1 = _mm_loadr_pd((double*)&(scratch[4].r));
            //src1 = _mm_load_pd((double*)&(scratch[4].r));
            //src1 = _mm_shuffle_pd(src1, src1, _MM_SHUFFLE2(0, 1));
            src2 = _mm_load_pd((double*)&(scratch[5].r));
#  if 1
            minus1 = _mm_set_pd(-0.0, 0.0);
            //minus2 = _mm_set_pd(0.0, -0.0);
            src3 = _mm_xor_pd(src1, minus1);
            //src4 = _mm_xor_pd(src1, minus2);
#  else
            minus1 = _mm_set_pd(-1.0, 1.0);
            //minus2 = _mm_set_pd(1.0, -1.0);
            src3 = _mm_mul_pd(src1, minus1);
            //src4 = _mm_mul_pd(src1, minus2);
#  endif
            src5 = _mm_add_pd(src2, src3);
            //src6 = _mm_add_pd(src2, src4);
            src6 = _mm_addsub_pd(src2, src1);  // a0-b0, a1+b1
            _mm_store_pd((double*)&(Fout[m].r), src5);
            _mm_store_pd((double*)&(Fout[m3].r), src6);
            }
# else
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
# endif
            ++Fout;
        }while(--k);
    }
# else
    do {
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);

        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );

        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
        ++Fout;
    }while(--k);
#endif
}

static void kf_bfly3(
         kiss_fft_cpx * Fout,
         const size_t fstride,
         const kiss_fft_cfg st,
         size_t m
         )
{
     size_t k=m;
     const size_t m2 = 2*m;
     kiss_fft_cpx *tw1,*tw2;
     __declspec(align(16)) kiss_fft_cpx scratch[5];
     kiss_fft_cpx epi3;
     epi3 = st->twiddles[fstride*m];

     tw1=tw2=st->twiddles;

     do{
         C_FIXDIV(*Fout,3); C_FIXDIV(Fout[m],3); C_FIXDIV(Fout[m2],3);

         C_MUL(scratch[1],Fout[m] , *tw1);
         C_MUL(scratch[2],Fout[m2] , *tw2);

         C_ADD(scratch[3],scratch[1],scratch[2]);
         C_SUB(scratch[0],scratch[1],scratch[2]);
         tw1 += fstride;
         tw2 += fstride*2;
#if 1
        {__m128d src1 = _mm_load_pd((double*)&((*Fout).r));
        __m128d src2 = _mm_load_pd((double*)&(scratch[3].r));
        __m128d src3 = _mm_set1_pd(-0.5);
        __m128d src4 = _mm_mul_pd(src2, src3);
        __m128d src5 = _mm_add_pd(src4, src1);
        _mm_store_pd((double*)&(Fout[m].r), src5);}
#else
         Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
         Fout[m].i = Fout->i - HALF_OF(scratch[3].i);
#endif
         C_MULBYSCALAR( scratch[0] , epi3.i );

         C_ADDTO(*Fout,scratch[3]);
#if 1
        {__m128d src1 = _mm_loadr_pd((double*)&(scratch[0].r));
        __m128d src2 = _mm_load_pd((double*)&(Fout[m].r));
# if 1
        __m128d minus1 = _mm_set_pd(-0.0, 0.0);
            //minus2 = _mm_set_pd(0.0, -0.0);
        __m128d src3 = _mm_xor_pd(src1, minus1);
            //src4 = _mm_xor_pd(src1, minus2);
# else
        __m128d minus1 = _mm_set_pd(-1.0, 1.0);
            //minus2 = _mm_set_pd(1.0, -1.0);
        __m128d src3 = _mm_mul_pd(src1, minus1);
            //src4 = _mm_mul_pd(src1, minus2);
# endif
        __m128d src5 = _mm_add_pd(src2, src3);
        __m128d src6 = _mm_addsub_pd(src2, src1);  // a0-b0, a1+b1
        _mm_store_pd((double*)&(Fout[m2].r), src5);
        _mm_store_pd((double*)&(Fout[m].r), src6);
        }
#else
         Fout[m2].r = Fout[m].r + scratch[0].i;
         Fout[m2].i = Fout[m].i - scratch[0].r;

         Fout[m].r -= scratch[0].i;
         Fout[m].i += scratch[0].r;
#endif
         ++Fout;
     }while(--k);
}
#if (1) // 0:original, 1:modified
static void kf_bfly5(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    int u;
    __declspec(align(16)) kiss_fft_cpx scratch[13];
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx *tw;
    kiss_fft_cpx ya,yb;
    ya = twiddles[fstride*m];
    yb = twiddles[fstride*2*m];

    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;

    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

        //C_ADD( scratch[7],scratch[1],scratch[4]);
        //C_SUB( scratch[10],scratch[1],scratch[4]);
        C_ADDSUB(scratch[7] ,scratch[10] ,scratch[1] ,scratch[4]);

        //C_ADD( scratch[8],scratch[2],scratch[3]);
        //C_SUB( scratch[9],scratch[2],scratch[3]);
        C_ADDSUB(scratch[8] ,scratch[9] ,scratch[2] ,scratch[3]);

        //Fout0->r += scratch[7].r + scratch[8].r;
        //Fout0->i += scratch[7].i + scratch[8].i;
        C_ADD2(*Fout0, scratch[7], scratch[8]);

        // scratch[5] = scratch[0] + S_MUL(scratch[7],ya.r) + S_MUL(scratch[8],yb.r);
#if 1
        {__m128d src1 = _mm_load_pd((double*)&(scratch[0].r));
        __m128d src2 = _mm_load_pd((double*)&(scratch[7].r));
        __m128d src3 = _mm_load_pd((double*)&(scratch[8].r));
        __m128d mag1 = _mm_load1_pd(&ya.r);
        __m128d mag2 = _mm_load1_pd(&yb.r);
        //__m128d mag1 = _mm_set1_pd(ya.r);
        //__m128d mag2 = _mm_set1_pd(yb.r);
# if (0)
        src1 = _mm_fmadd_pd(src1, src2, mag1);
        src1 = _mm_fmadd_pd(src1, src3, mag2);
# else
        src2 = _mm_mul_pd(src2, mag1);  // S_MUL(scratch[7], ya.r)
        src3 = _mm_mul_pd(src3, mag2);  // S_MUL(scratch[8], yb.r)
        src1 = _mm_add_pd(src1, src2);  // scratch[0] += S_MUL(scratch[7], ya.r)
        src1 = _mm_add_pd(src1, src3);  // scratch[0] += S_MUL(scratch[8], yb.r)
# endif
        _mm_store_pd((double*)&(scratch[5].r), src1);}
#else
        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r,ya.r) + S_MUL(scratch[8].r,yb.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i,ya.r) + S_MUL(scratch[8].i,yb.r);
#endif
#if 1
        // temp1 = 10 x ya.i
        // temp2 = 9 x yb.i
        // 6 = temp1 + temp2
        // 6 = (6.i, 6.r)
        // 6.i = -6.i;
        {__m128d src1 = _mm_load_pd((double*)&(scratch[9].r));
        __m128d src2 = _mm_load_pd((double*)&(scratch[10].r));
        __m128d mag1 = _mm_load1_pd(&ya.i);
        __m128d mag2 = _mm_load1_pd(&yb.i);
        //__m128d mag1 = _mm_set1_pd(ya.i);
        //__m128d mag2 = _mm_set1_pd(yb.i);
        __m128d src3 = _mm_mul_pd(src1, mag2);  // S_MUL(scratch[9], yb.i)
        __m128d src4 = _mm_mul_pd(src2, mag1);  // S_MUL(scratch[10], ya.i)
        __m128d src5 = _mm_add_pd(src3, src4);  // S_MUL(scratch[9], yb.i) + S_MUL(scratch[10], ya.i)
        __m128d minus1 = _mm_set_pd(0.0, -0.0);
        src1 = _mm_xor_pd(src5, minus1);
        src1 = _mm_shuffle_pd(src1, src1, _MM_SHUFFLE2(0, 1));
        _mm_store_pd((double*)&(scratch[6].r), src1);}
#else
        scratch[6].r =  S_MUL(scratch[10].i,ya.i) + S_MUL(scratch[9].i,yb.i);
        scratch[6].i = -S_MUL(scratch[10].r,ya.i) - S_MUL(scratch[9].r,yb.i);
#endif
        //C_SUB(*Fout1,scratch[5],scratch[6]);
        //C_ADD(*Fout4,scratch[5],scratch[6]);
        C_ADDSUB(*Fout4 ,*Fout1 ,scratch[5] ,scratch[6]);
#if 1
        {__m128d src1 = _mm_load_pd((double*)&(scratch[0].r));
        __m128d src2 = _mm_load_pd((double*)&(scratch[7].r));
        __m128d src3 = _mm_load_pd((double*)&(scratch[8].r));
        __m128d mag1 = _mm_load1_pd(&yb.r);
        __m128d mag2 = _mm_load1_pd(&ya.r);
        //__m128d mag1 = _mm_set1_pd(yb.r);
        //__m128d mag2 = _mm_set1_pd(ya.r);
        src2 = _mm_mul_pd(src2, mag1);  // S_MUL(scratch[7], yb.r)
        src3 = _mm_mul_pd(src3, mag2);  // S_MUL(scratch[8], ya.r)
        src1 = _mm_add_pd(src1, src2);  // scratch[0] += S_MUL(scratch[7], yb.r)
        src1 = _mm_add_pd(src1, src3);  // scratch[0] += S_MUL(scratch[8], ya.r)
        _mm_store_pd((double*)&(scratch[11].r), src1);}
#else
        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r,yb.r) + S_MUL(scratch[8].r,ya.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i,yb.r) + S_MUL(scratch[8].i,ya.r);
#endif
#if 1
        {__m128d src1 = _mm_load_pd((double*)&(scratch[9].r));
        __m128d src2 = _mm_load_pd((double*)&(scratch[10].r));
        __m128d mag1 = _mm_load1_pd(&ya.i);
        __m128d mag2 = _mm_load1_pd(&yb.i);
        //__m128d mag1 = _mm_set1_pd(ya.i);
        //__m128d mag2 = _mm_set1_pd(yb.i);
        __m128d src3 = _mm_mul_pd(src1, mag1);  // S_MUL(scratch[9], ya.i)
        __m128d src4 = _mm_mul_pd(src2, mag2);  // S_MUL(scratch[10], yb.i)
        __m128d src5 = _mm_sub_pd(src3, src4);  // S_MUL(scratch[9], yb.i) - S_MUL(scratch[10], yb.i)
        __m128d minus1 = _mm_set_pd(0.0, -0.0);
        src1 = _mm_xor_pd(src5, minus1);
        src1 = _mm_shuffle_pd(src1, src1, _MM_SHUFFLE2(0, 1));
        _mm_store_pd((double*)&(scratch[12].r), src1);}
#else
        scratch[12].r = - S_MUL(scratch[10].i,yb.i) + S_MUL(scratch[9].i,ya.i);
        scratch[12].i = S_MUL(scratch[10].r,yb.i) - S_MUL(scratch[9].r,ya.i);
#endif
        //C_ADD(*Fout2,scratch[11],scratch[12]);
        //C_SUB(*Fout3,scratch[11],scratch[12]);
        C_ADDSUB(*Fout2 ,*Fout3 ,scratch[11] ,scratch[12]);

        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
}
#else
static void kf_bfly5(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    int u;
    __declspec(align(16)) kiss_fft_cpx scratch[13];
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx *tw;
    kiss_fft_cpx ya,yb;
    ya = twiddles[fstride*m];
    yb = twiddles[fstride*2*m];

    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;

    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

        C_ADD( scratch[7],scratch[1],scratch[4]);
        C_SUB( scratch[10],scratch[1],scratch[4]);

        C_ADD( scratch[8],scratch[2],scratch[3]);
        C_SUB( scratch[9],scratch[2],scratch[3]);

        Fout0->r += scratch[7].r + scratch[8].r;
        Fout0->i += scratch[7].i + scratch[8].i;

        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r,ya.r) + S_MUL(scratch[8].r,yb.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i,ya.r) + S_MUL(scratch[8].i,yb.r);

        scratch[6].r =  S_MUL(scratch[10].i,ya.i) + S_MUL(scratch[9].i,yb.i);
        scratch[6].i = -S_MUL(scratch[10].r,ya.i) - S_MUL(scratch[9].r,yb.i);

        C_SUB(*Fout1,scratch[5],scratch[6]);
        C_ADD(*Fout4,scratch[5],scratch[6]);

        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r,yb.r) + S_MUL(scratch[8].r,ya.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i,yb.r) + S_MUL(scratch[8].i,ya.r);

        scratch[12].r = - S_MUL(scratch[10].i,yb.i) + S_MUL(scratch[9].i,ya.i);
        scratch[12].i = S_MUL(scratch[10].r,yb.i) - S_MUL(scratch[9].r,ya.i);

        C_ADD(*Fout2,scratch[11],scratch[12]);
        C_SUB(*Fout3,scratch[11],scratch[12]);

        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
}
#endif
/* perform the butterfly for one stage of a mixed radix FFT */
static void kf_bfly_generic(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m,
        int p
        )
{
    int u,k,q1,q;
    kiss_fft_cpx * twiddles = st->twiddles;
    __declspec(align(16)) kiss_fft_cpx t;
    int Norig = st->nfft;

    kiss_fft_cpx * scratch = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC(sizeof(kiss_fft_cpx)*p);

    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratch[q1] = Fout[ k  ];
            C_FIXDIV(scratch[q1],p);
            k += m;
        }

        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratch[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig) twidx-=Norig;
                C_MUL(t,scratch[q] , twiddles[twidx] );
                C_ADDTO( Fout[ k ] ,t);
            }
            k += m;
        }
    }
    KISS_FFT_TMP_FREE(scratch);
}

static
void kf_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        const size_t fstride,
        int in_stride,
        int * factors,
        const kiss_fft_cfg st
        )
{
    kiss_fft_cpx * Fout_beg=Fout;
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */
    const kiss_fft_cpx * Fout_end = Fout + p*m;

#ifdef _OPENMP
    // use openmp extensions at the 
    // top-level (not recursive)
    if (fstride==1 && p<=5)
    {
        int k;

        // execute the p different work units in different threads
#       pragma omp parallel for
        for (k=0;k<p;++k) 
            kf_work( Fout +k*m, f+ fstride*in_stride*k,fstride*p,in_stride,factors,st);
        // all threads have joined by this point

        switch (p) {
            case 2: kf_bfly2(Fout,fstride,st,m); break;
            case 3: kf_bfly3(Fout,fstride,st,m); break; 
            case 4: kf_bfly4(Fout,fstride,st,m); break;
            case 5: kf_bfly5(Fout,fstride,st,m); break; 
            default: kf_bfly_generic(Fout,fstride,st,m,p); break;
        }
        return;
    }
#endif

    if (m==1) {
        do{
            *Fout = *f;
            f += fstride*in_stride;
        }while(++Fout != Fout_end );
    }else{
        do{
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m, 
            // each one takes a decimated version of the input
            kf_work( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
        }while( (Fout += m) != Fout_end );
    }

    Fout=Fout_beg;

    // recombine the p smaller DFTs 
    switch (p) {
        case 2: kf_bfly2(Fout,fstride,st,m); break;
        case 3: kf_bfly3(Fout,fstride,st,m); break; 
        case 4: kf_bfly4(Fout,fstride,st,m); break;
        case 5: kf_bfly5(Fout,fstride,st,m); break; 
        default: kf_bfly_generic(Fout,fstride,st,m,p); break;
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where 
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static 
void kf_factor(int n,int * facbuf)
{
    int p=4;
    double floor_sqrt;
    floor_sqrt = floor( sqrt((double)n) );

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p > floor_sqrt)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
kiss_fft_cfg kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem )
{
    kiss_fft_cfg st=NULL;
    size_t memneeded = sizeof(struct kiss_fft_state)
        + sizeof(kiss_fft_cpx)*(nfft-1); /* twiddle factors*/

    if ( lenmem==NULL ) {
        st = ( kiss_fft_cfg)KISS_FFT_MALLOC( memneeded );
    }else{
        if (mem != NULL && *lenmem >= memneeded)
            st = (kiss_fft_cfg)mem;
        *lenmem = memneeded;
    }
    if (st) {
        int i;
        st->nfft=nfft;
        st->inverse = inverse_fft;

        for (i=0;i<nfft;++i) {
            const double pi=3.141592653589793238462643383279502884197169399375105820974944;
            double phase = -2*pi*i / nfft;
            if (st->inverse)
                phase *= -1;
            kf_cexp(st->twiddles+i, phase );
        }

        kf_factor(nfft,st->factors);
    }
    return st;
}


void kiss_fft_stride(kiss_fft_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int in_stride)
{
    if (fin == fout) {
        //NOTE: this is not really an in-place FFT algorithm.
        //It just performs an out-of-place FFT into a temp buffer
        kiss_fft_cpx * tmpbuf = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC( sizeof(kiss_fft_cpx)*st->nfft);
        kf_work(tmpbuf,fin,1,in_stride, st->factors,st);
        memcpy(fout,tmpbuf,sizeof(kiss_fft_cpx)*st->nfft);
        KISS_FFT_TMP_FREE(tmpbuf);
    }else{
        kf_work( fout, fin, 1,in_stride, st->factors,st );
    }
}

void kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    kiss_fft_stride(cfg,fin,fout,1);
}


void kiss_fft_cleanup(void)
{
    // nothing needed any more
}

int kiss_fft_next_fast_size(int n)
{
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}
