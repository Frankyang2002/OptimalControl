function ode_stance = get_dyn_stance(t,in2,in3)
%GET_DYN_STANCE
%    ODE_STANCE = GET_DYN_STANCE(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    23-Nov-2022 15:51:46

dq1 = in2(5,:);
dq2 = in2(6,:);
dy = in2(4,:);
q1 = in2(2,:);
q2 = in2(3,:);
u1 = in3(1,:);
u2 = in3(2,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = sin(q1);
t5 = sin(q2);
t6 = q1+q2;
t7 = dq1.^2;
t8 = dq2.^2;
t9 = q1.*2.0;
t10 = q1.*3.0;
t11 = q2.*2.0;
t12 = q1.*4.0;
t13 = q2.*3.0;
t14 = q1.*5.0;
t15 = q2.*4.0;
t16 = q1.*6.0;
t17 = q2.*5.0;
t36 = -q2;
t40 = -u2;
t18 = cos(t9);
t19 = cos(t10);
t20 = cos(t11);
t21 = cos(t12);
t22 = cos(t13);
t23 = cos(t14);
t24 = cos(t15);
t25 = t2.^2;
t26 = t3.^2;
t27 = sin(t9);
t28 = sin(t10);
t29 = sin(t11);
t30 = sin(t12);
t31 = sin(t13);
t32 = sin(t14);
t33 = sin(t15);
t34 = cos(t6);
t35 = sin(t6);
t37 = -t11;
t38 = -t13;
t39 = -t15;
t41 = q2+t6;
t42 = q1+t6;
t43 = t6+t11;
t44 = t6+t9;
t45 = t6+t13;
t46 = t6+t10;
t47 = t6+t15;
t48 = t6+t12;
t66 = t2.*2.11e+2;
t67 = t4.*2.11e+2;
t68 = q1+t36;
t72 = t6.*2.0;
t76 = t6.*3.0;
t82 = t6.*4.0;
t87 = t6.*5.0;
t98 = t2.*(4.0./2.5e+1);
t100 = t3.*6.58e+2;
t101 = t4.*(4.0./2.5e+1);
t102 = t9+t36;
t104 = t10+t36;
t107 = t2.*t3.*4.22e+2;
t158 = dq1.*dq2.*t5.*3.968e-4;
t159 = dq1.*dq2.*t5.*7.936e-4;
t162 = t5.*t7.*3.968e-4;
t164 = t2.*1.655928e-1;
t173 = t2.*7.353843043764298e+23;
t174 = t3.*9.552528378623778e+25;
t186 = t3.*1.421426683960568e+49;
t195 = dq2.*dy.*t4.*4.666436787997162e+50;
t198 = dq2.*dy.*t2.*2.891112236337437e+50;
t199 = t5.*1.240168945141381e+53;
t202 = dq1.*dy.*t4.*2.78431323224793e+51;
t204 = t3.*3.524159239043114e+53;
t209 = dq1.*dy.*t2.*1.999949107735983e+51;
t213 = dq1.*dq2.*t4.*1.323058483100174e+52;
t221 = dq1.*dq2.*t2.*1.383449314884857e+52;
t230 = t4.*t8.*7.574799303084676e+51;
t256 = t2.*t7.*1.983133712640686e+52;
t261 = t4.*t7.*1.801636956303656e+52;
t263 = t2.*u1.*2.177162285889733e+54;
t264 = t2.*u2.*2.188758905295069e+54;
t266 = t2.*t8.*7.793787542442714e+51;
t297 = t4.*u1.*8.647501485276129e+54;
t325 = t4.*u2.*7.285528714784229e+54;
t49 = cos(t41);
t50 = cos(t42);
t51 = cos(t43);
t52 = cos(t44);
t53 = cos(t45);
t54 = cos(t46);
t55 = cos(t47);
t56 = cos(t48);
t57 = t34.^2;
t58 = sin(t41);
t59 = sin(t42);
t60 = sin(t43);
t61 = sin(t44);
t62 = sin(t45);
t63 = sin(t46);
t64 = sin(t47);
t65 = sin(t48);
t69 = q1+t37;
t70 = q1+t38;
t71 = q1+t39;
t73 = t6+t41;
t74 = t6+t42;
t75 = t6+t43;
t77 = t6+t44;
t78 = t41+t72;
t79 = t42+t72;
t80 = t6+t46;
t81 = t43+t72;
t83 = t44+t72;
t84 = t6+t48;
t86 = t46+t72;
t88 = t44+t76;
t89 = cos(t68);
t93 = sin(t68);
t97 = t35.*6.2e+1;
t99 = -t66;
t103 = t9+t37;
t105 = t10+t37;
t106 = t34.*2.73e+2;
t108 = cos(t72);
t112 = cos(t76);
t118 = cos(t82);
t123 = cos(t87);
t125 = sin(t72);
t129 = sin(t76);
t135 = sin(t82);
t140 = sin(t87);
t142 = t34.*(4.0./2.5e+1);
t143 = -t107;
t144 = t35.*(4.0./2.5e+1);
t145 = cos(t102);
t147 = cos(t104);
t150 = sin(t102);
t152 = sin(t104);
t154 = t3.*t34.*6.2e+1;
t156 = t34.*t66;
t157 = t2.*t34.*-2.11e+2;
t161 = dq2.*dy.*t35.*4.96e-3;
t167 = t34.*2.43288e-2;
t168 = -t164;
t172 = t34.*7.258760166127491e+22;
t175 = t3.*t34.*7.090928421933952e+22;
t176 = -t174;
t177 = t26.*4.66583090163254e+25;
t179 = t25.*1.551660882234267e+26;
t181 = t2.*t34.*3.063196790105801e+25;
t184 = t2.*t3.*t34.*2.992371794056128e+25;
t187 = t21.*2.03008617538409e+48;
t188 = t22.*5.439487853197403e+48;
t189 = t24.*1.045436988055579e+48;
t190 = t20.*3.387990991939753e+49;
t192 = t18.*2.727530954127463e+49;
t200 = dq1.*dq2.*t23.*2.886323673154288e+49;
t203 = dq1.*dq2.*t32.*1.443161836577144e+49;
t206 = t24.*8.255351938159233e+50;
t207 = t33.*2.323027642837929e+51;
t208 = dq2.*dy.*t28.*8.883791172695913e+49;
t212 = dq2.*dy.*t35.*5.364238157337525e+50;
t214 = t27.*9.023474847613485e+52;
t215 = dq1.*dy.*t23.*2.318919588010254e+50;
t217 = dq2.*dy.*t34.*2.941960469588733e+50;
t218 = -t204;
t219 = dq2.*dy.*t19.*1.117431868474281e+50;
t223 = t22.*1.016658396699888e+52;
t225 = t31.*1.011886810454755e+52;
t233 = -t213;
t250 = t7.*t23.*2.607551955325683e+50;
t255 = dq1.*dy.*t19.*2.054473029750693e+51;
t258 = dq1.*dq2.*t28.*1.550836958489879e+51;
t259 = t18.*2.908966132934357e+53;
t260 = -t221;
t262 = dq1.*dy.*t34.*1.954465651117035e+51;
t265 = dq1.*dy.*t35.*2.124517625842464e+51;
t267 = dq1.*dq2.*t19.*2.124482335677761e+51;
t270 = t19.*u1.*3.777199760521161e+53;
t271 = t19.*u2.*3.893165954574528e+53;
t272 = -t230;
t280 = t8.*t19.*1.260524701596498e+51;
t289 = t7.*t19.*5.03939613834198e+51;
t305 = t20.*2.029084792075241e+52;
t310 = -t256;
t312 = dq1.*dy.*t32.*1.159459794005127e+50;
t315 = t34.*u1.*1.780516673754214e+53;
t316 = -t261;
t323 = -t264;
t324 = -t266;
t327 = dq1.*dq2.*t35.*2.316501239229616e+52;
t330 = dq1.*dq2.*t34.*4.237894003305549e+52;
t341 = t8.*t28.*9.214808112123509e+50;
t347 = t35.*u1.*1.120597049408335e+54;
t348 = t28.*u2.*9.559141944808468e+53;
t364 = dq1.*dy.*t28.*1.209627489934543e+51;
t371 = -t297;
t378 = t7.*t35.*1.87362985550639e+52;
t381 = t29.*1.667259794858647e+53;
t383 = t7.*t32.*1.303775977662842e+50;
t388 = t34.*u2.*4.22784486267587e+54;
t390 = t35.*u2.*4.63016366383503e+54;
t397 = t7.*t34.*3.242548324750888e+52;
t408 = t28.*u1.*1.098881834214756e+54;
t409 = t8.*t34.*2.55496268591265e+52;
t417 = t8.*t35.*1.321738797552873e+52;
t424 = t30.*6.882433062429596e+51;
t452 = t21.*1.942876945730539e+52;
t599 = t7.*t28.*3.111447913648438e+51;
t85 = t6.*2.0+t74;
t90 = cos(t69);
t91 = cos(t70);
t92 = cos(t71);
t94 = sin(t69);
t95 = sin(t70);
t96 = sin(t71);
t109 = cos(t73);
t110 = cos(t74);
t111 = cos(t75);
t113 = cos(t77);
t114 = cos(t78);
t115 = cos(t79);
t116 = cos(t80);
t117 = cos(t81);
t119 = cos(t83);
t120 = cos(t84);
t122 = cos(t86);
t124 = cos(t88);
t126 = sin(t73);
t127 = sin(t74);
t128 = sin(t75);
t130 = sin(t77);
t131 = sin(t78);
t132 = sin(t79);
t133 = sin(t80);
t134 = sin(t81);
t136 = sin(t83);
t137 = sin(t84);
t139 = sin(t86);
t141 = sin(t88);
t146 = cos(t103);
t148 = cos(t105);
t149 = t57.*3.1e+1;
t151 = sin(t103);
t153 = sin(t105);
t160 = t67+t97;
t163 = -t161;
t165 = t98+t142;
t166 = t101+t144;
t169 = -t167;
t178 = -t175;
t182 = t57.*1.192933294743937e+25;
t183 = t99+t106+t143+t154;
t185 = -t184;
t191 = -t188;
t193 = -t190;
t194 = t54.*1.100121739622618e+48;
t196 = t118.*9.569421518428202e+45;
t197 = -t192;
t210 = t50.*9.6226165597707e+48;
t224 = -t206;
t226 = -t207;
t227 = -t208;
t231 = dq1.*dy.*t123.*6.552818316031e+47;
t232 = dq2.*dy.*t123.*6.552818316031e+47;
t234 = dq1.*dq2.*t53.*7.149582608244554e+49;
t235 = dq1.*dq2.*t55.*5.198417882217304e+49;
t236 = dq1.*dy.*t53.*5.922459917282918e+48;
t237 = dq1.*dy.*t55.*1.30504214616334e+49;
t238 = dq2.*dy.*t53.*6.604936303452963e+49;
t239 = dq2.*dy.*t55.*1.30504214616334e+49;
t240 = t7.*t64.*3.26260536540835e+48;
t241 = dq1.*dq2.*t62.*2.264429487462966e+49;
t242 = dq1.*dq2.*t64.*2.599208941108652e+49;
t243 = dq1.*dy.*t62.*7.562768538679051e+48;
t244 = dq1.*dy.*t64.*6.525210730816699e+48;
t245 = dq2.*dy.*t62.*2.842314293722722e+49;
t246 = dq2.*dy.*t64.*6.525210730816699e+48;
t247 = t108.*7.455341961461138e+48;
t248 = -t214;
t249 = dq2.*dy.*t65.*3.406943470536393e+49;
t251 = dq2.*dy.*t56.*6.813886941072785e+49;
t252 = dq1.*dq2.*t140.*6.397874072379621e+47;
t253 = dq1.*dy.*t140.*3.2764091580155e+47;
t254 = dq2.*dy.*t140.*3.2764091580155e+47;
t257 = -t219;
t269 = t145.*1.820602245999661e+48;
t273 = dq1.*dq2.*t123.*1.279574814475924e+48;
t276 = dq2.*dy.*t89.*7.976547596497567e+49;
t278 = t7.*t55.*6.525210730816699e+48;
t279 = t8.*t55.*3.240854662972294e+49;
t291 = t8.*t64.*1.620427331486147e+49;
t293 = dq2.*dy.*t58.*2.621912753379068e+50;
t303 = dq1.*dq2.*t65.*3.819002927050442e+50;
t304 = t7.*t123.*3.2764091580155e+47;
t306 = -t250;
t309 = -t255;
t311 = dq1.*dy.*t65.*1.286712966878907e+50;
t314 = dq1.*dq2.*t60.*6.771858007986333e+50;
t317 = dq1.*dy.*t56.*2.573425933757814e+50;
t318 = dq2.*dy.*t49.*1.374432887537354e+50;
t319 = dq2.*dy.*t60.*2.200356687619272e+50;
t322 = dq2.*dy.*t61.*3.141201821572785e+50;
t328 = t8.*t123.*2.966520670712743e+47;
t329 = dq1.*dy.*t89.*6.733036430226725e+49;
t332 = dq2.*dy.*t52.*5.40426821955609e+50;
t333 = t52.*u1.*4.295695129363725e+51;
t334 = dq1.*dy.*t60.*3.800533472411481e+50;
t335 = t63.*5.700423799737386e+52;
t336 = dq2.*dy.*t51.*4.148474847920439e+50;
t342 = t7.*t60.*5.590203322001241e+50;
t343 = dq1.*dy.*t51.*7.691995153046382e+50;
t345 = -t270;
t349 = t8.*t60.*6.764153160759472e+50;
t356 = t7.*t53.*3.508109869322277e+50;
t357 = t8.*t53.*2.483346465936417e+50;
t359 = t8.*t56.*4.283481830060431e+50;
t362 = dq1.*dq2.*t56.*7.638005854100884e+50;
t363 = dq1.*dy.*t58.*1.537546905791802e+51;
t366 = t7.*t62.*1.715049524595959e+50;
t367 = t8.*t62.*1.241673232968208e+50;
t374 = dq1.*dq2.*t51.*1.274692213695904e+51;
t379 = -t305;
t380 = dq2.*dy.*t147.*4.99199150896981e+48;
t384 = dq1.*dq2.*t61.*5.916306526156541e+51;
t385 = t7.*t56.*5.359019672931933e+50;
t387 = -t315;
t391 = dq1.*dy.*t61.*1.167807044726878e+51;
t396 = dq1.*dy.*t52.*2.051720736602323e+51;
t414 = -t327;
t418 = -t330;
t420 = dq1.*dy.*t49.*2.578036825754562e+51;
t422 = t56.*u2.*7.274001192425121e+52;
t425 = t50.*3.071638092209698e+53;
t429 = dq2.*dy.*t93.*5.193466174507533e+50;
t431 = t7.*t51.*1.27391968560338e+51;
t433 = dq1.*dq2.*t52.*1.145491872659541e+52;
t436 = -t347;
t438 = dq1.*dq2.*t129.*3.743588297737053e+50;
t439 = -t348;
t440 = t8.*t51.*1.34282262860623e+51;
t451 = -t364;
t459 = -t378;
t462 = dq1.*dy.*t93.*1.556089038770322e+51;
t464 = dq1.*dy.*t112.*4.426416434034435e+50;
t465 = t7.*t140.*1.63820457900775e+47;
t466 = dq1.*dq2.*t49.*1.245675360950397e+52;
t467 = dq1.*dq2.*t58.*6.830084138667289e+51;
t469 = -t383;
t472 = dq1.*dq2.*t112.*7.352235964905973e+50;
t474 = t61.*u1.*1.690357748287585e+53;
t477 = -t390;
t482 = -t397;
t490 = dq2.*dy.*t112.*2.935407651272702e+50;
t492 = dq2.*dy.*t152.*5.523047770860085e+49;
t493 = -t409;
t495 = t8.*t58.*3.09580387888305e+51;
t496 = -t417;
t498 = t8.*t61.*3.369131615185034e+51;
t502 = t8.*t140.*1.483260335356371e+47;
t503 = dq2.*dy.*t129.*1.503368895679472e+50;
t504 = dq1.*dy.*t129.*2.233890667739242e+50;
t505 = t51.*u1.*4.965908931107528e+53;
t506 = t53.*u1.*9.811934818684385e+52;
t507 = t51.*u2.*3.726394913403458e+53;
t508 = t53.*u2.*9.811934818684385e+52;
t509 = t54.*3.686618759989555e+52;
t512 = t7.*t61.*4.456429502539522e+51;
t514 = t62.*u1.*1.962386963736877e+53;
t515 = t62.*u2.*1.962386963736877e+53;
t516 = t65.*u2.*1.454800238485024e+53;
t520 = t8.*t49.*6.275473583537027e+51;
t531 = -t452;
t532 = dq1.*dq2.*t93.*8.160843908625319e+51;
t540 = t8.*t65.*2.141740915030215e+50;
t542 = t7.*t58.*7.898099007993095e+51;
t543 = t7.*t49.*1.59373815221467e+52;
t544 = dq1.*dq2.*t147.*4.111110458200527e+50;
t550 = dq1.*dy.*t152.*2.629640895065318e+50;
t552 = t108.*3.247319028412848e+53;
t559 = t8.*t112.*2.349239932644793e+50;
t561 = t7.*t112.*5.671725861531781e+50;
t563 = t7.*t65.*2.679509836465966e+50;
t564 = t7.*t89.*2.170043331392813e+51;
t567 = t8.*t52.*6.627464760515255e+51;
t568 = dq1.*dq2.*t152.*5.654156981375906e+50;
t575 = t8.*t93.*3.585466464384815e+51;
t578 = t7.*t52.*8.453792126195312e+51;
t587 = dq1.*dq2.*t89.*6.894112336759805e+51;
t589 = t89.*u1.*5.081875125160895e+53;
t591 = t49.*u1.*3.087645687011528e+54;
t593 = t49.*u2.*3.762288247497702e+54;
t595 = t52.*u2.*1.013351161758973e+54;
t598 = t125.*6.699427126061682e+53;
t600 = t145.*4.326173777486791e+52;
t603 = dq1.*dy.*t147.*1.860915817002883e+50;
t605 = t60.*u1.*1.016375025032179e+54;
t609 = t60.*u2.*8.309854301258076e+53;
t610 = t61.*u2.*1.72693756653906e+54;
t613 = t7.*t147.*1.668267817880145e+50;
t616 = t7.*t93.*7.644306747348541e+51;
t620 = t112.*u1.*1.780516673754214e+53;
t623 = t112.*u2.*1.506462947826021e+53;
t625 = t58.*u1.*6.174236405409992e+54;
t631 = t129.*u2.*3.065020236180744e+53;
t641 = t93.*u1.*1.109435491232995e+54;
t642 = t8.*t147.*5.184003669107339e+50;
t644 = t129.*u1.*3.583223245792043e+53;
t652 = t89.*u2.*4.072125241087548e+54;
t660 = t93.*u2.*4.228257920443376e+54;
t662 = t58.*u2.*7.024189481653535e+54;
t665 = t59.*4.145431683168164e+53;
t666 = t8.*t89.*6.147690560476166e+51;
t671 = t8.*t129.*1.179267068482598e+50;
t672 = t135.*4.239497445330838e+51;
t674 = t147.*u2.*4.504599879763673e+53;
t677 = t150.*8.733080202928303e+51;
t678 = t152.*u2.*6.099599282557298e+53;
t681 = t7.*t129.*2.798932747083513e+50;
t693 = t118.*2.101560242803112e+51;
t700 = t8.*t152.*1.691479995506761e+50;
t703 = t7.*t152.*6.324935785900953e+50;
t121 = cos(t85);
t138 = sin(t85);
t155 = -t149;
t170 = (dq1.*dy.*t160)./1.25e+4;
t201 = t113.*3.428854090497467e+47;
t205 = t115.*1.108635950376181e+47;
t211 = t111.*2.000422676302168e+47;
t220 = t109.*1.679181674777088e+48;
t228 = -t210;
t268 = t146.*5.083273555670504e+48;
t274 = -t231;
t275 = -t232;
t277 = dq2.*dy.*t90.*9.901583228965647e+48;
t281 = dq1.*dq2.*t91.*4.717511489225447e+49;
t282 = dq1.*dq2.*t92.*2.588333589890624e+49;
t283 = -t234;
t284 = dq1.*dy.*t92.*4.441353110330076e+49;
t285 = -t236;
t286 = dq2.*dy.*t91.*6.008043945806426e+49;
t287 = -t237;
t288 = -t239;
t290 = t7.*t96.*2.73988687138056e+48;
t294 = dq1.*dq2.*t95.*6.65104418708844e+49;
t295 = dq1.*dq2.*t96.*1.946687868026982e+49;
t296 = -t241;
t298 = dq1.*dy.*t96.*2.220676555165038e+49;
t299 = -t243;
t300 = dq2.*dy.*t95.*1.463583227242023e+49;
t301 = -t244;
t302 = -t246;
t307 = -t253;
t308 = -t254;
t320 = dq1.*dy.*t134.*4.048848565421163e+48;
t321 = dq2.*dy.*t134.*4.048848565421163e+48;
t326 = t124.*1.111266682617775e+50;
t331 = t141.*2.222533365235551e+50;
t338 = -t269;
t340 = dq2.*dy.*t94.*8.946335114325431e+49;
t350 = dq1.*dy.*t119.*5.724590861012852e+49;
t351 = dq1.*dq2.*t131.*2.749896478428539e+49;
t352 = dq1.*dq2.*t119.*8.35236687958643e+49;
t353 = dq2.*dy.*t133.*2.779770352093204e+49;
t355 = t7.*t92.*7.0296867002207e+49;
t358 = -t279;
t360 = dq1.*dy.*t91.*2.848613362368667e+50;
t365 = t7.*t95.*9.260730834926507e+49;
t368 = -t291;
t369 = t8.*t119.*5.622495377375961e+48;
t372 = dq1.*dy.*t95.*1.113662058168548e+50;
t375 = -t303;
t376 = dq1.*dy.*t114.*8.284670668526472e+48;
t382 = dq1.*dy.*t136.*2.862295430506426e+49;
t386 = -t314;
t389 = t7.*t117.*4.048848565421163e+48;
t392 = dq1.*dq2.*t117.*1.834785599280021e+49;
t393 = dq1.*dy.*t117.*8.097697130842326e+48;
t394 = dq2.*dy.*t117.*8.097697130842326e+48;
t395 = dq2.*dy.*t114.*3.807303455684881e+49;
t398 = t7.*t134.*2.024424282710581e+48;
t399 = t8.*t134.*3.100655148268361e+48;
t400 = dq2.*dy.*t119.*3.771616406774711e+49;
t401 = dq1.*dq2.*t114.*5.691197644855405e+49;
t402 = dq2.*dy.*t116.*5.559540704186407e+49;
t403 = dq1.*dq2.*t134.*9.173927996400105e+48;
t404 = dq2.*dy.*t131.*1.85962696964775e+49;
t405 = -t322;
t406 = dq2.*dy.*t136.*1.885808203387355e+49;
t407 = dq1.*dy.*t133.*9.197601495879203e+49;
t410 = dq1.*dq2.*t90.*7.410398945661114e+50;
t411 = t120.*5.137772522233711e+51;
t412 = t122.*1.122505960214959e+51;
t415 = -t328;
t416 = -t329;
t419 = t139.*2.245011920429917e+51;
t421 = -t332;
t426 = dq1.*dq2.*t136.*4.176183439793215e+49;
t427 = t8.*t136.*2.811247688687981e+48;
t430 = -t342;
t432 = t8.*t114.*7.304016171515008e+49;
t435 = t7.*t119.*5.576601743768712e+49;
t437 = dq2.*dy.*t127.*2.259179721194898e+50;
t441 = dq1.*dq2.*t116.*3.566351459055939e+50;
t443 = t7.*t91.*3.947444058650678e+50;
t445 = t8.*t91.*2.386580914647564e+50;
t446 = dq1.*dq2.*t94.*1.580104601752489e+51;
t447 = -t359;
t449 = t8.*t94.*2.755396351370328e+50;
t450 = -t362;
t453 = t8.*t95.*1.208631700625067e+50;
t456 = t8.*t133.*8.651645685332139e+49;
t457 = -t374;
t463 = t8.*t131.*3.652008085757504e+49;
t468 = dq2.*dy.*t110.*3.942200757530361e+50;
t470 = dq1.*dy.*t94.*9.24089370977988e+50;
t471 = -t385;
t473 = dq1.*dy.*t90.*1.698160374193243e+51;
t476 = t8.*t117.*6.201310296536721e+48;
t478 = -t391;
t481 = -t396;
t484 = t132.*3.018041571146393e+52;
t485 = dq1.*dq2.*t153.*8.937126161451127e+49;
t487 = dq1.*dy.*t153.*7.201580764004825e+49;
t489 = dq1.*dy.*t116.*1.839520299175841e+50;
t491 = t109.*4.236859954281581e+52;
t494 = t111.*9.74653423703833e+51;
t497 = t7.*t114.*6.519008224563678e+49;
t499 = t128.*2.025350414918948e+52;
t500 = t137.*1.027554504446742e+52;
t510 = t7.*t136.*2.788300871884356e+49;
t511 = dq1.*dq2.*t133.*1.783175729527969e+50;
t518 = -t431;
t524 = dq1.*dy.*t127.*8.028417018165593e+50;
t525 = t7.*t116.*3.323222623926985e+50;
t526 = dq1.*dq2.*t110.*4.304391507976758e+51;
t533 = t114.*u1.*2.583205358356322e+52;
t534 = t116.*u1.*7.274001192425121e+52;
t535 = t119.*u1.*1.58923145347004e+52;
t536 = t114.*u2.*2.583205358356322e+52;
t537 = t116.*u2.*8.86323264589516e+52;
t538 = t119.*u2.*1.431899552544429e+52;
t545 = t131.*u1.*5.166410716712643e+52;
t546 = t136.*u1.*3.17846290694008e+52;
t547 = t131.*u2.*5.166410716712643e+52;
t549 = t136.*u2.*2.863799105088858e+52;
t551 = dq1.*dy.*t131.*3.702087752316331e+48;
t554 = -t464;
t555 = -t466;
t556 = dq1.*dy.*t110.*1.523792277795244e+51;
t557 = t7.*t90.*4.360114543499975e+51;
t558 = -t467;
t565 = dq1.*dq2.*t148.*1.210160497659368e+50;
t566 = dq1.*dy.*t148.*2.488808331299388e+50;
t569 = t7.*t153.*1.735545397446302e+49;
t571 = t8.*t90.*1.023016802401118e+51;
t572 = -t490;
t573 = dq1.*dq2.*t127.*2.218969552150403e+51;
t574 = -t492;
t576 = t146.*8.306540072465183e+51;
t579 = t7.*t94.*1.103531989199209e+51;
t580 = -t495;
t582 = t8.*t116.*1.730329137066428e+50;
t583 = t151.*1.235665505769331e+52;
t585 = -t502;
t586 = -t503;
t588 = -t504;
t590 = t90.*u1.*3.523406397895234e+53;
t592 = -t505;
t594 = t91.*u2.*9.811934818684385e+52;
t596 = -t508;
t597 = -t509;
t602 = t7.*t133.*1.661611311963493e+50;
t606 = t94.*u2.*4.255584911929738e+53;
t607 = -t514;
t608 = t95.*u2.*1.304603623936806e+53;
t611 = -t516;
t612 = -t520;
t617 = t8.*t127.*1.058377638164763e+51;
t618 = t7.*t127.*2.305217984728003e+51;
t619 = -t532;
t626 = -t540;
t627 = -t542;
t628 = -t543;
t629 = -t544;
t630 = t133.*u1.*1.454800238485024e+53;
t633 = t133.*u2.*1.772646529179032e+53;
t635 = -t550;
t638 = -t559;
t639 = -t563;
t640 = t7.*t148.*3.698968828958756e+50;
t645 = t8.*t110.*2.115916395053157e+51;
t646 = -t575;
t650 = t115.*1.480112496166049e+52;
t651 = -t591;
t653 = t90.*u2.*8.489315329002763e+53;
t654 = -t595;
t655 = -t598;
t657 = -t603;
t658 = t7.*t131.*3.251851284672057e+49;
t659 = t94.*u1.*1.109023420284983e+54;
t664 = -t609;
t667 = t7.*t110.*4.579433472370616e+51;
t668 = t130.*1.674830011679983e+53;
t669 = t126.*9.36937489159519e+52;
t670 = -t616;
t673 = t110.*u1.*9.347421356961021e+53;
t675 = t110.*u2.*1.128686117606224e+54;
t676 = -t623;
t679 = t127.*u2.*2.209960450028199e+54;
t683 = -t641;
t684 = t113.*8.361902467269735e+52;
t685 = -t642;
t686 = t127.*u1.*1.85449836982885e+54;
t687 = -t644;
t689 = -t652;
t692 = -t662;
t694 = -t665;
t696 = -t671;
t697 = -t672;
t699 = -t678;
t704 = t172+t173+t178;
t705 = t176+t177+t179+t181+t182+t185-3.369145534876864e+26;
t171 = -t170;
t180 = t100+t155+t157+3.29e+2;
t216 = -t201;
t222 = -t205;
t229 = -t211;
t292 = dq2.*dy.*t121.*5.738597922405914e+48;
t313 = dq2.*dy.*t138.*2.869298961202957e+48;
t337 = dq1.*dq2.*t138.*4.432049532314928e+48;
t339 = dq1.*dy.*t121.*3.508525882595364e+48;
t344 = t7.*t121.*1.737114110188021e+48;
t346 = dq1.*dq2.*t121.*8.864099064629856e+48;
t354 = -t277;
t361 = -t286;
t370 = -t294;
t373 = -t298;
t377 = t8.*t121.*4.714909587302337e+48;
t413 = -t326;
t428 = t8.*t138.*2.357454793651169e+48;
t442 = -t352;
t444 = -t355;
t448 = -t360;
t454 = t121.*u1.*1.57331900925611e+51;
t455 = t121.*u2.*1.57331900925611e+51;
t458 = -t376;
t460 = t138.*u1.*3.146638018512219e+51;
t461 = t138.*u2.*3.146638018512219e+51;
t475 = -t389;
t479 = -t392;
t480 = -t395;
t483 = -t398;
t486 = -t403;
t488 = -t404;
t501 = -t419;
t513 = -t426;
t517 = dq1.*dy.*t138.*1.754262941297682e+48;
t519 = -t432;
t521 = -t435;
t522 = t7.*t138.*8.685570550940105e+47;
t523 = -t437;
t527 = -t441;
t528 = -t445;
t529 = -t446;
t530 = -t449;
t541 = -t456;
t553 = -t463;
t560 = -t468;
t562 = -t473;
t570 = -t487;
t577 = -t494;
t581 = -t497;
t584 = -t500;
t601 = -t510;
t604 = -t511;
t614 = -t524;
t615 = -t525;
t621 = -t533;
t622 = -t534;
t624 = -t535;
t632 = -t547;
t634 = -t549;
t636 = -t551;
t637 = -t556;
t643 = -t569;
t647 = -t576;
t648 = -t579;
t649 = -t582;
t656 = -t602;
t661 = -t606;
t663 = -t608;
t680 = -t633;
t682 = -t640;
t688 = -t650;
t690 = -t653;
t691 = -t658;
t695 = -t669;
t698 = -t675;
t701 = -t684;
t702 = -t686;
t706 = 1.0./t705;
t423 = -t337;
t434 = -t346;
t539 = -t455;
t548 = -t460;
t707 = t186+t187+t189+t191+t193+t194+t196+t197+t216+t220+t222+t228+t229+t247+t268+t338+4.607443903526238e+49;
t708 = 1.0./t707;
t709 = t198+t199+t200+t209+t215+t217+t225+t226+t235+t238+t248+t251+t257+t260+t262+t267+t273+t274+t275+t276+t278+t280+t281+t282+t283+t284+t285+t287+t288+t289+t292+t304+t306+t309+t310+t317+t318+t324+t325+t331+t335+t336+t339+t343+t344+t350+t354+t356+t357+t358+t361+t369+t371+t377+t380+t381+t393+t394+t400+t401+t402+t408+t410+t415+t416+t418+t420+t421+t424+t433+t434+t436+t439+t440+t442+t443+t444+t447+t448+t450+t457+t458+t461+t471+t472+t474+t475+t476+t477+t479+t480+t481+t482+t484+t489+t493+t499+t501+t515+t518+t519+t521+t526+t527+t528+t545+t546+t548+t554+t555+t557+t560+t561+t562+t564+t565+t566+t567+t571+t572+t578+t581+t583+t584+t587+t605+t607+t610+t611+t612+t613+t615+t625+t628+t629+t630+t631+t632+t634+t637+t638+t645+t649+t655+t657+t659+t660+t661+t663+t664+t666+t667+t668+t677+t679+t680+t682+t683+t685+t687+t692+t694+t695+t697+t699+t702;
t710 = t195+t202+t203+t212+t218+t223+t224+t227+t233+t240+t242+t245+t249+t252+t258+t259+t263+t265+t271+t272+t290+t293+t295+t296+t299+t300+t301+t302+t307+t308+t311+t312+t313+t316+t319+t320+t321+t323+t333+t334+t340+t341+t345+t349+t351+t353+t363+t365+t366+t367+t368+t370+t372+t373+t375+t379+t382+t384+t386+t387+t388+t399+t405+t406+t407+t411+t412+t413+t414+t422+t423+t425+t427+t428+t429+t430+t438+t451+t453+t454+t459+t462+t465+t469+t470+t478+t483+t485+t486+t488+t491+t496+t498+t506+t507+t512+t513+t517+t522+t523+t529+t530+t531+t536+t537+t538+t539+t541+t552+t553+t558+t568+t570+t573+t574+t577+t580+t585+t586+t588+t589+t590+t592+t593+t594+t596+t597+t599+t600+t601+t604+t614+t617+t618+t619+t620+t621+t622+t624+t626+t627+t635+t636+t639+t643+t646+t647+t648+t651+t654+t656+t670+t673+t674+t676+t681+t688+t689+t690+t691+t693+t696+t698+t700+t701+t703-4.805394706774209e+53;
t711 = (t35.*t708.*t709)./1.5625e+5;
t713 = (t166.*t708.*t709)./2.5e+4;
t714 = (t708.*t710)./1.25e+4;
t715 = (t34.*t708.*t710)./7.8125e+4;
t717 = t165.*t708.*t710.*(-8.0e-5);
t712 = -t711;
t716 = t165.*t714;
t718 = t163+t171+t714+3.22749;
t720 = t159+t162+t168+t169+t713+t717+u1;
t719 = t40+t158+t167+t712+t715;
ode_stance = [dy;dq1;dq2;t706.*(t3.*7.258760166127491e+22-t26.*3.545464210966976e+22+2.560140983948985e+23).*(t161+t170-t714-3.22749).*-4.0e+3-t183.*t706.*t719.*1.42962266571249e+25+t704.*t706.*t720.*1.25e+4;t706.*t720.*(t57.*3.545464210966976e+22-1.146641877440026e+24).*1.5625e+5+t704.*t706.*(t161+t170-t714-3.22749).*1.25e+4-t180.*t706.*t719.*1.787028332140613e+26;t183.*t706.*(t161+t170-t714-3.22749).*1.42962266571249e+25-t706.*t719.*(t3.*-4.0796e+4+t25.*4.4521e+4+t57.*9.61e+2+t2.*t34.*1.3082e+4-1.00016e+5).*5.764607523034235e+24+t180.*t706.*t720.*1.787028332140613e+26];