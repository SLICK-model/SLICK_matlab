function [Y_mc] = SLICK_simulation(K_y, G, y0, dt, Nf, M_n, nt, k_max)

%%

%  Simulating the SLICK model by solving the stochastic differentiation equaitons


%       Inputs:  
%
%          K_y:  Inflated Koopman operator
%            G:  De-whitening filter
%           y0:  Initial condition
%           dt:  Time step
%           Nf:  Frequency bins
%          M_n:  Modal rank
%           nt:  Size of the training set
%        k_max:  Total number of Monte-Carlo simulations


%      Outputs:  
%
%         Y_mc:  Inflated state vectors for all the simulations with size k_max*(2*Nf*M_n)*(nt+1)





%%

w1           =       zeros(2*(Nf)*M_n,1);
g            =       [zeros((Nf)*M_n,(Nf)*2*M_n);zeros((Nf)*M_n) G ];
Y_mc         =       zeros(k_max,2*(Nf)*M_n,nt+1);
seed         =       roundn(1e5*rand +1e4*rand+1e3*rand+1e2*rand,0);


for k = 1:k_max
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['k = ' num2str(k)])
    
    Y0       =   zeros((Nf)*M_n*2,nt+1);
    
    Y0(:,1)  =   y0;
    
    
    for it = 1:nt
        
        if mod(it,100)==0
            disp([' it= ' num2str(it) ', t=' num2str(it*dt)])
        end
        
        rng(seed);
        
        for j=1:1:2*(Nf)*M_n
            [ n1, seed ] = r8_normal_01 ( seed );
            w1(j) = n1 * sqrt ( 1 / (dt));
        end
          
         Y0(:,it+1)=Y0(:,it)+(K_y*Y0(:,it)  + 1* g*w1)*(dt);
           
        
    end
    
    Y_mc(k,:,:)         =    Y0;
    
    seed  = roundn(1e5*rand +1e4*rand+1e3*rand+1e2*rand,0);
    
end


end





function [ x, seed ] = r8_normal_01 ( seed )

%*****************************************************************************80
%
%% R8_NORMAL_01 returns a unit pseudonormal R8.
%
%  Discussion:
%
%    The standard normal probability distribution function (PDF) has
%    mean 0 and standard deviation 1.
%
%    Because this routine uses the Box Muller method, it requires pairs
%    of uniform random values to generate a pair of normal random values.
%    This means that on every other call, the code can use the second
%    value that it calculated.
%
%    However, if the user has changed the SEED value between calls,
%    the routine automatically resets itself and discards the saved data.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    13 August 2008
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer SEED, a seed for the random number generator.
%
%    Output, real X, a sample of the standard normal PDF.
%
%    Output, integer SEED, an updated seed for the random number generator.
%
  persistent seed1;
  persistent seed2;
  persistent seed3;
  persistent used;
  persistent v2;

  if ( size ( used ) == 0 )
    used = 0;
    seed1 = 0;
    seed2 = 0;
    seed3 = 0;
    v2 = 0.0;
  end
%
%  If USED is odd, but the input SEED does not match
%  the output SEED on the previous call, then the user has changed
%  the seed.  Wipe out internal memory.
%
  if ( mod ( used, 2 ) == 1 )

    if ( seed ~= seed2 )
      used = 0;
      seed1 = 0;
      seed2 = 0;
      seed3 = 0;
      v2 = 0.0;
    end

  end
%
%  If USED is even, generate two uniforms, create two normals,
%  return the first normal and its corresponding seed.
%
  if ( mod ( used, 2 ) == 0 )

    seed1 = seed;

    [ r1, seed ] = r8_uniform_01 ( seed );

    if ( r1 == 0.0 )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'R8_NORMAL_01 - Fatal error!\n' );
      fprintf ( 1, '  R8_UNIFORM_01 returned a value of 0.\n' );
      error ( 'R8_NORMAL_01 - Fatal error!\n' );
    end

    seed2 = seed;

    [ r2, seed ] = r8_uniform_01 ( seed );

    seed3 = seed;

    v1 = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );
    v2 = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * pi * r2 );

    x = v1;
    seed = seed2;
%
%  If USED is odd (and the input SEED matched the output value from
%  the previous call), return the second normal and its corresponding seed.
%
  else

    x = v2;
    seed = seed3;

  end

  used = used + 1;

  return
end


function [ r, seed ] = r8_uniform_01 ( seed )

%*****************************************************************************80
%
%% R8_UNIFORM_01 returns a unit pseudorandom R8.
%
%  Discussion:
%
%    This routine implements the recursion
%
%      seed = 16807 * seed mod ( 2**31 - 1 )
%      r8_uniform_01 = seed / ( 2**31 - 1 )
%
%    The integer arithmetic never requires more than 32 bits,
%    including a sign bit.
%
%    If the initial seed is 12345, then the first three computations are
%
%      Input     Output      R8_UNIFORM_01
%      SEED      SEED
%
%         12345   207482415  0.096616
%     207482415  1790989824  0.833995
%    1790989824  2035175616  0.947702
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license. 
%
%  Modified:
%
%    21 September 2006
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Paul Bratley, Bennett Fox, Linus Schrage,
%    A Guide to Simulation,
%    Second Edition,
%    Springer, 1987,
%    ISBN: 0387964673,
%    LC: QA76.9.C65.B73.
%
%    Bennett Fox,
%    Algorithm 647:
%    Implementation and Relative Efficiency of Quasirandom
%    Sequence Generators,
%    ACM Transactions on Mathematical Software,
%    Volume 12, Number 4, December 1986, pages 362-376.
%
%    Pierre L'Ecuyer,
%    Random Number Generation,
%    in Handbook of Simulation,
%    edited by Jerry Banks,
%    Wiley, 1998,
%    ISBN: 0471134031,
%    LC: T57.62.H37.
%
%    Peter Lewis, Allen Goodman, James Miller,
%    A Pseudo-Random Number Generator for the System/360,
%    IBM Systems Journal,
%    Volume 8, Number 2, 1969, pages 136-143.
%
%  Parameters:
%
%    Input, integer SEED, the integer "seed" used to generate
%    the output random number.  SEED should not be 0.
%
%    Output, real R, a random value between 0 and 1.
%
%    Output, integer SEED, the updated seed.  This would
%    normally be used as the input seed on the next call.
%
  i4_huge = 2147483647;

  if ( seed == 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8_UNIFORM_01 - Fatal error!\n' );
    fprintf ( 1, '  Input SEED = 0!\n' );
    error ( 'R8_UNIFORM_01 - Fatal error!' );
  end

  seed = floor ( seed );

  seed = mod ( seed, i4_huge );

  if ( seed < 0 ) 
    seed = seed + i4_huge;
  end 

  k = floor ( seed / 127773 );

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
    seed = seed + i4_huge;
  end

  r = seed * 4.656612875E-10;

  return
end



