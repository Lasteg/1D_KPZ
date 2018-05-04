c ============================          
         function ranpang(iseed)
c Pang, p.47         
c ============================          
         i2e30 = 2**30
         ia=16807
         ic=i2e30-1    ! ic=2**31-1, but 2**31 is a overflow
         ic=ic+i2e30
          
         iq=ic/ia
         ir=mod(ic,ia)

         ih=iseed/iq
         il=mod(iseed,iq)

         it=ia*il-ir*ih
         if(it.gt.0) then
           iseed=it
         else
           iseed=it+ic
         endif

         ranpang=iseed/float(ic)

         return
         end
