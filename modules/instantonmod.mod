  �,  c   k820309    [          18.0        U5�[                                                                                                          
       instantonmod.f90 INSTANTONMOD              PARTITION                                                     
                        �                                  
                                                                %         @                                                    
       #X                                                                  
               &                   &                                                                                                                                                                                                                         
                &                                                                                        	     
                                                    
            #         @                                                      #X    #GRAD                                                                  
               &                   &                                                                                                         
               &                   &                                           #         @                                                      #X    #HESS                                                                  
               &                   &                                                                                                         
               &                   &                   &                   &                                                                                                    #         @ �                                                    #A              
D@                                                  
 -              &                                           %         @                                                     
       #X    #A    #B              D @                                                  
               &                   &                   &                                                     D @                                                  
               &                   &                                                     D @                                                  
               &                   &                                           #         @                                                       #X    #A    #B    #ANSWER              D @                                                  
               &                   &                   &                                                                                                         
               &                   &                                                                                                         
               &                   &                                                     D                                                    
               &                   &                   &                                           #         @                                                       #X    #ANSWER              D @                                                  
 	              &                   &                   &                                                     D                                                    
 
              &                   &                                           #         @                                                        #ATOMSIN !   #THETA1 "   #THETA2 #   #THETA3 $   #ORIGIN %                                              !                   
               &                   &                                                     D @                               "     
                 D @                               #     
                 D                                 $     
                D                                 %                    
     p          5 r        5 r                      #         @                                  &                    #ATOMS '   #AXIS (   #THETA )             D                                 '                   
               &                   &                                                                                      (                       @                               )     
       #         @                                   *                    #ATOMSIN +   #THETA1 ,   #THETA2 -   #THETA3 .   #ORIGIN /   #ATOMSOUT 0                                              +                   
               &                   &                                                     D @                               ,     
                 D @                               -     
                 D @                               .     
                                                 /                    
     p          5 r        5 r                                D @                               0                   
               &                   &                                           #         @                                   1                    #VEC 2   #AXIS 3   #THETA 4             D                                 2                   
               &                                                                                      3                       @                               4     
       %         @                                5                           #N1 6   #N2 7   #N3 8   #STRING 9             
                                  6                     
                                  7                     
                                  8                     
                                9                    1 #         @                                   :                    #X ;   #Y <   #YP1 =   #YPN >   #Y2 ?          0  
 @                               ;                   
              &                                                     
 @                               <                   
              &                                                     
                                  =     
                
                                  >     
                D@                               ?                   
               &                                           #         @                                  @                    #A A   #B B   #C C   #R D   #U E             
 @                               A                   
              &                                                  0  
 @                               B                   
               &                                                     
 @                               C                   
 !             &                                                     
 @                               D                   
 "             &                                                     D@                               E                   
 #              &                                           %         @                                F                           #NN G   #STRING H             
                                  G                    %             &                                                     
                                H                    1 %         @                                 I                    
       #XA J   #YA K   #Y2A L   #X M             
 @                               J                   
 &             &                                                     
 @                               K                   
 '             &                                                     
 @                               L                   
 (             &                                                     
  @                               M     
      %         @                                N                           #XX O   #X P             
 @                               O                   
 ,             &                                                     
                                  P     
      %         @                                Q                    
       #XA R   #YA S   #Y2A T   #X U             
 @                               R                   
 )             &                                                     
 @                               S                   
 *             &                                                     
 @                               T                   
 +             &                                                     
  @                               U     
      #         @                                   V                    #X W   #COM X                                              W                   
 /              &                   &                                                     D                                 X                   
 0              &                                           %         @                                 Y                    
       #X1 Z   #X2 [   #LAMPATH \   #PATH ]   #SPLINEPATH ^              @                               Z     
                  @                               [     
                  @                               \                   
 1              &                                                      @                               ]                   
 2              &                                                      @                               ^                   
 3              &                                           %         @                                 _                    
       #X1 `   #X2 a                                              `                   
 4              &                   &                                                                                      a                   
 5              &                   &                                              �   &      fn#fn "   �      b   uapp(INSTANTONMOD    �   @   J   POTENTIAL       @   j   GENERAL    `  @       N+GENERAL    �  W       POT+POTENTIAL     �  �   a   POT%X+POTENTIAL    �  @       NDIM+GENERAL    �  @       NATOM+GENERAL      �       MASS+GENERAL    �  @       BETAN+GENERAL "   �  @       FIXEDENDS+GENERAL !   '  Y       VPRIME+POTENTIAL #   �  �   a   VPRIME%X+POTENTIAL &   $  �   a   VPRIME%GRAD+POTENTIAL '   �  Y       VDOUBLEPRIME+POTENTIAL )   !  �   a   VDOUBLEPRIME%X+POTENTIAL ,   �  �   a   VDOUBLEPRIME%HESS+POTENTIAL    �  @       NDOF+GENERAL    �  O       QSORTC    (  �   a   QSORTC%A    �  e       UM    	  �   a   UM%X    �	  �   a   UM%A    y
  �   a   UM%B      i       UMPRIME    �  �   a   UMPRIME%X    B  �   a   UMPRIME%A    �  �   a   UMPRIME%B    �  �   a   UMPRIME%ANSWER    F  [       UMHESSIAN    �  �   a   UMHESSIAN%X !   ]  �   a   UMHESSIAN%ANSWER      �       GET_ALIGN "   �  �   a   GET_ALIGN%ATOMSIN !   *  @   a   GET_ALIGN%THETA1 !   j  @   a   GET_ALIGN%THETA2 !   �  @   a   GET_ALIGN%THETA3 !   �  �   a   GET_ALIGN%ORIGIN    ~  h       ROTATE_ATOMS #   �  �   a   ROTATE_ATOMS%ATOMS "   �  @   a   ROTATE_ATOMS%AXIS #   �  @   a   ROTATE_ATOMS%THETA    
  �       ALIGN_ATOMS $   �  �   a   ALIGN_ATOMS%ATOMSIN #   A  @   a   ALIGN_ATOMS%THETA1 #   �  @   a   ALIGN_ATOMS%THETA2 #   �  @   a   ALIGN_ATOMS%THETA3 #     �   a   ALIGN_ATOMS%ORIGIN %   �  �   a   ALIGN_ATOMS%ATOMSOUT    9  f       ROTATE_VEC    �  �   a   ROTATE_VEC%VEC     +  @   a   ROTATE_VEC%AXIS !   k  @   a   ROTATE_VEC%THETA    �  t       ASSERT_EQ      @   a   ASSERT_EQ%N1    _  @   a   ASSERT_EQ%N2    �  @   a   ASSERT_EQ%N3 !   �  L   a   ASSERT_EQ%STRING    +  p       SPLINE    �  �   a   SPLINE%X    '  �   a   SPLINE%Y    �  @   a   SPLINE%YP1    �  @   a   SPLINE%YPN    3  �   a   SPLINE%Y2    �  k       TRIDAG    *  �   a   TRIDAG%A    �  �   a   TRIDAG%B    B  �   a   TRIDAG%C    �  �   a   TRIDAG%R    Z  �   a   TRIDAG%U    �  d       ASSERT_EQN    J   �   a   ASSERT_EQN%NN "   �   L   a   ASSERT_EQN%STRING    "!  p       SPLINT    �!  �   a   SPLINT%XA    "  �   a   SPLINT%YA    �"  �   a   SPLINT%Y2A    6#  @   a   SPLINT%X    v#  _       LOCATE    �#  �   a   LOCATE%XX    a$  @   a   LOCATE%X    �$  p       SPLIN_GRAD    %  �   a   SPLIN_GRAD%XA    �%  �   a   SPLIN_GRAD%YA    )&  �   a   SPLIN_GRAD%Y2A    �&  @   a   SPLIN_GRAD%X    �&  X       CENTREOFMASS    M'  �   a   CENTREOFMASS%X !   �'  �   a   CENTREOFMASS%COM    }(  �       FINDMIDDLE    )  @   a   FINDMIDDLE%X1    D)  @   a   FINDMIDDLE%X2 #   �)  �   a   FINDMIDDLE%LAMPATH     *  �   a   FINDMIDDLE%PATH &   �*  �   a   FINDMIDDLE%SPLINEPATH    (+  `       EUCLIDDIST    �+  �   a   EUCLIDDIST%X1    ,,  �   a   EUCLIDDIST%X2 