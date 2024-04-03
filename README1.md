<div align="center">
  	<h1>
    	BdryReach User Guide
  	</h1>
	<br />
    <div>
		<a href="https://github.com/BdryReach/BdryReach">
			<img width="428" src="result_picture/BdryReach.svg" alt="BdryReach System">
		</a>
	</div>
</div>

## 1. Artifact Requirements

To run BdrtReach in a Linux system, it is necessary to install the **cmake** tool and the various **third-party libraries** required for BdryReach.

### 1.1 Resource requirements

Our implementation utilizes the floating point linear programming solver GLPK and C++ library Eigen. We adopt the approach outlined in  \cite{althoff2008reachability} to compute outer-approximations appeared in our method. All  experiments herein are run on Ubuntu 20.04.3 LTS in virtual machine 
with CPU 12th Gen Intel Core i9-12900K × 8  and RAM 15.6 GB. 
| Library | Website | Version |
| --- | --- | --- |
| git | [https://git-scm.com/](https://git-scm.com/) | 2.25.1 |
| Cmake | [https://cmake.org/](https://cmake.org/) | latest |
| Eigen | [http://eigen.tuxfamily.org/index.php?title=Main Page](http://eigen.tuxfamily.org/index.php?title=Main%20Page) | 3.34 |
| Python | [https://www.python.org/](https://www.python.org/) | 2.7.18 |
| Capd | [http://capd.ii.uj.edu.pl/](http://capd.ii.uj.edu.pl/) | latest |
| boost | [https://www.boost.org/](https://www.boost.org/) | 1.67.0.0 |
| GLPK | [https://www.gnu.org/software/glpk/](https://www.gnu.org/software/glpk/) | 4.65-2 |

### 1.2 Evaluation runtime

**We provided all the experiments in our paper for the artifact evaluation. The examples files are located in ./examples, which contains all the examples showcased in Table 3-5. It has tow folders, ./examples/BdryReach and ./examples/CORA, which include examples run by BdryReach and CORA respectively.
Here are the computation times of each file, which are tested in the environment stated in 1.1.**

<body lang=ZH-CN style='tab-interval:21.0pt;word-wrap:break-word'>

<div class=WordSection1 style='layout-grid:15.6pt'>

<div align=center>

<table class=MsoTableGrid border=1 cellspacing=0 cellpadding=0
 style='border-collapse:collapse;border:none;mso-border-alt:solid windowtext .5pt;
 mso-yfti-tbllook:1184;mso-padding-alt:0cm 5.4pt 0cm 5.4pt'>
 <tr style='mso-yfti-irow:0;mso-yfti-firstrow:yes;height:14.45pt'>
  <td width=282 colspan=2 style='width:211.25pt;border:solid windowtext 1.0pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=GramE><span
  lang=EN-US>..</span></span><span lang=EN-US>/Table3_examples</span></p>
  </td>
  <td width=286 colspan=2 style='width:214.15pt;border:solid windowtext 1.0pt;
  border-left:none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=GramE><span
  lang=EN-US>..</span></span><span lang=EN-US>/Table4_examples</span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:1;height:13.85pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US>Flie</span></span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Runtime</span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>File</span></p>
  </td>
  <td width=140 style='width:104.95pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Runtime</span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:2;height:14.45pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>ElectroOsc_T3</span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>23.56</span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>ElectroOsc_T4</span></p>
  </td>
  <td width=140 style='width:104.95pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>73.76</span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:3;height:13.85pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Roessler_T3</span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>27.72</span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Roessler_T4</span></p>
  </td>
  <td width=140 style='width:104.95pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>63.42</span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:4;height:14.45pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>LotkaVolterra_T3</span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>10.43</span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>LotkaVolterra_T4</span></p>
  </td>
  <td width=140 style='width:104.95pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>81.81</span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:5;height:13.85pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Tank6_T3</span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>50.83</span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Tank6_T4</span></p>
  </td>
  <td width=140 style='width:104.95pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>129.58</span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:6;height:14.45pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>BiologicalSysteml_T3</span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>1.74</span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>BiologicalSysteml_T4</span></p>
  </td>
  <td width=140 style='width:104.95pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>462.87</span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:7;height:13.85pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>BiologicalSystemll_T3</span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>72.47</span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>BiologicalSystemll_T4</span></p>
  </td>
  <td width=140 style='width:104.95pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>261.78</span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:8;mso-yfti-lastrow:yes;height:14.45pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Tank12_T3</span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>235.88</span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Tank12_T4</span></p>
  </td>
  <td width=140 style='width:104.95pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>377.85</span></p>
  </td>
 </tr>
</table>

</div>

<p class=MsoNormal><span lang=EN-US><o:p>&nbsp;</o:p></span></p>

<div align=center>

<table class=MsoTableGrid border=1 cellspacing=0 cellpadding=0 width=652
 style='width:489.05pt;border-collapse:collapse;border:none;mso-border-alt:
 solid windowtext .5pt;mso-yfti-tbllook:1184;mso-padding-alt:0cm 5.4pt 0cm 5.4pt'>
 <tr style='mso-yfti-irow:0;mso-yfti-firstrow:yes;height:14.45pt'>
  <td width=319 colspan=2 style='width:239.4pt;border:solid windowtext 1.0pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=GramE><span
  lang=EN-US>..</span></span><span lang=EN-US>/Table5_examples_longTime<o:p></o:p></span></p>
  </td>
  <td width=333 colspan=2 style='width:249.65pt;border:solid windowtext 1.0pt;
  border-left:none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=GramE><span
  lang=EN-US>..</span></span><span lang=EN-US>/Table4_examples<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:1;height:13.85pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US>Flie</span></span><span lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 style='width:83.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Runtime<o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>File<o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Runtime<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:2;height:14.45pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>ElectroOsc_biginit_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 style='width:83.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>11.29</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>ElectroOsc_biginit_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>4.79<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:3;height:13.85pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>Roessler_biginit_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 style='width:83.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>24.01</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Roessler_biginit_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>15.88<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:4;height:14.45pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>LotkaVolterra_biginit_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 style='width:83.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>64.17</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>LotkaVolterra_biginit_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>15.45<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:5;height:13.85pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>Tank6_biginit_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 style='width:83.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>80.91</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Tank6_biginit_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>66.83<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:6;height:16.3pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:16.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:
  white'>BiologicalSysteml_biginit</span></span><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'> _T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 style='width:83.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:16.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>281.05</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:16.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>BiologicalSystems_biginit</span></span><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'> _T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:16.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>118.65<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:7;height:13.85pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:
  white'>BiologicalSystemll_biginit</span></span><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'> _T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 style='width:83.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>142.02</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>BiologicalSystemsl_biginit</span></span><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'> _T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>142.17<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:8;mso-yfti-lastrow:yes;height:20.55pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:20.55pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>Tank12_big_init_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 style='width:83.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:20.55pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>235.57</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:20.55pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Tank12_big_init_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:20.55pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>162.46<o:p></o:p></span></p>
  </td>
 </tr>
</table>

</div>

<p class=MsoNormal><span lang=EN-US><o:p>&nbsp;</o:p></span></p>

<div align=center>

<table class=MsoTableGrid border=1 cellspacing=0 cellpadding=0
 style='border-collapse:collapse;border:none;mso-border-alt:solid windowtext .5pt;
 mso-yfti-tbllook:1184;mso-padding-alt:0cm 5.4pt 0cm 5.4pt'>
 <tr style='mso-yfti-irow:0;mso-yfti-firstrow:yes;height:14.45pt'>
  <td width=282 colspan=2 style='width:211.25pt;border:solid windowtext 1.0pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=GramE><span
  lang=EN-US>..</span></span><span lang=EN-US>/Table3_examples<o:p></o:p></span></p>
  </td>
  <td width=286 colspan=2 style='width:214.15pt;border:solid windowtext 1.0pt;
  border-left:none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=GramE><span
  lang=EN-US>..</span></span><span lang=EN-US>/Table4_examples<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:1;height:13.85pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US>Flie</span></span><span lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Runtime<o:p></o:p></span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>File<o:p></o:p></span></p>
  </td>
  <td width=140 style='width:104.95pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Runtime<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:2;height:14.45pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>ElectroOsc_T3<o:p></o:p></span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>36.50</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>ElectroOsc_T4</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=140 rowspan=7 style='width:104.95pt;border-top:none;border-left:
  none;border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>a
  few seconds<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:3;height:13.85pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Roessler_T3<o:p></o:p></span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>36.63</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Roessler_T4</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:4;height:14.45pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>LotkaVolterra_T3<o:p></o:p></span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>335.06</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>LotkaVolterra_T4</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:5;height:13.85pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Tank6_T3<o:p></o:p></span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>201.05</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Tank6_T4</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:6;height:14.45pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>BiologicalSysteml_T3<o:p></o:p></span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>125.73</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>BiologicalSysteml_T4</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:7;height:13.85pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>BiologicalSystemll_T3<o:p></o:p></span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>188.25</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>BiologicalSystemll_T4</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:8;mso-yfti-lastrow:yes;height:14.45pt'>
  <td width=142 style='width:106.35pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Tank12_T3<o:p></o:p></span></p>
  </td>
  <td width=140 style='width:104.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>1834.65</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=146 style='width:109.2pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Tank12_T4</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
</table>

</div>

<p class=MsoNormal><span lang=EN-US><o:p>&nbsp;</o:p></span></p>

<div align=center>

<table class=MsoTableGrid border=1 cellspacing=0 cellpadding=0 width=652
 style='width:489.05pt;border-collapse:collapse;border:none;mso-border-alt:
 solid windowtext .5pt;mso-yfti-tbllook:1184;mso-padding-alt:0cm 5.4pt 0cm 5.4pt'>
 <tr style='mso-yfti-irow:0;mso-yfti-firstrow:yes;height:14.45pt'>
  <td width=319 colspan=2 style='width:239.4pt;border:solid windowtext 1.0pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=GramE><span
  lang=EN-US>..</span></span><span lang=EN-US>/Table5_examples_longTime<o:p></o:p></span></p>
  </td>
  <td width=333 colspan=2 style='width:249.65pt;border:solid windowtext 1.0pt;
  border-left:none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=GramE><span
  lang=EN-US>..</span></span><span lang=EN-US>/Table4_examples<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:1;height:13.85pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US>Flie</span></span><span lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 style='width:83.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Runtime<o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>File<o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>Runtime<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:2;height:14.45pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>ElectroOsc_biginit_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=111 rowspan=7 style='width:83.55pt;border-top:none;border-left:
  none;border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US>a
  few seconds<o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>ElectroOsc_biginit_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>24.32</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:3;height:13.85pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>Roessler_biginit_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Roessler_biginit_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>36.54</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:4;height:14.45pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>LotkaVolterra_biginit_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>LotkaVolterra_biginit_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:14.45pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>153.57</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:5;height:13.85pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>Tank6_biginit_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Tank6_biginit_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>463.28</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:6;height:16.3pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:16.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:
  white'>BiologicalSysteml_biginit</span></span><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'> _T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:16.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>BiologicalSystems_biginit</span></span><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'> _T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:16.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>615.89</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:7;height:13.85pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:
  white'>BiologicalSystemll_biginit</span></span><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'> _T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span class=SpellE><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>BiologicalSystemsl_biginit</span></span><span
  lang=EN-US style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'> _T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:13.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>1494.38</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:8;mso-yfti-lastrow:yes;height:20.55pt'>
  <td width=208 style='width:155.85pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:20.55pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D;background:white'>Tank12_big_init_T5l</span></span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=210 style='width:157.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:20.55pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>Tank12_big_init_T5s</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
  <td width=123 style='width:92.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:20.55pt'>
  <p class=MsoNormal align=center style='text-align:center'><span lang=EN-US
  style='font-family:"Segoe UI",sans-serif;color:#0D0D0D'>1693.68</span><span
  lang=EN-US><o:p></o:p></span></p>
  </td>
 </tr>
</table>

## 2. Structure and Content

### 2.1 Outer-approximation and Inner-approximation of Reachable Set Computation Interface

### 2.1.1  Outer-approximation of Reachable Set Computation Interface 
```cpp
template <typename Number>
static vector<ReachableSet<Number>> BdReach(NonlinearSys<Number> mysys, ReachOptions<Number> options, Zonotope<Number> R0)
```
**Parameters:**
* **mysys:** differential equation for computing reachable sets.
* **options:** configuration for outer-approximation of reachable set computation.
* **R0:** initial set.



### 2.2 Test Case for Outer-approximation of Reachable Set Computation
**As an example, we perform the outer-approximation of the reachable set computation for the VanderPol model. The file computes the outer-approximation from the initial region ([1.23, 1.57], [2.34, 2.46]) over the time interval 0 - 6.74 seconds.The specific file location is:**
```RobotFramework
/examples/overVanderPol.cpp.
```
### 2.2.1 Include Files
```cpp
#include <overApprox/overApprox.h> // Header File with Interfaces for Computing Reachable Set Outer-approximation.
#include <plotter/matplotlibcpp.h> // Header file for Matplotlib C++ plotting library
#include <plotter/plotter.h> // Header file for result plotting
```
### 2.2.2 Definition of Differential Equations
**We define the form of differential equations using the Capd library. For detailed information on the differential equation system in Capd, please refer to the [Capd documentation](https://capd.sourceforge.net/capdDynSys/docs/html/maps.html) on ordinary differential equation systems.**


```cpp
double mu = 1.;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = in[1];
    out[1] = mu * (1-in[0]*in[0])*in[1] - in[0]/*+ in[2]*/;
}

// Input dimension of the differential equation
int dimIn = 3; // The input dimension of the differential equation. Since the default input includes control u, the input dimension is one greater than the output dimension.

// Output dimension of the differential equation
int dimOut = 2; // The output dimension of the differential equation.

// Parameter settings for the differential equation. Since this differential equation has no parameters, it is set to 0.
int noParam = 0;

// Maximum order for Taylor expansion of the differential equation
int MaxDerivativeOrder = 3; // The maximum order to which the differential equation is expanded using Taylor series.

// Creating IMap for interval computations
IMap f(_f, dimIn, dimOut, noParam, MaxDerivativeOrder); // Constructing IMap for interval Computations


```
### 2.2.3 Parameter Configuration for Computing Reachable Sets
**Here, we adopt the same parameter definitions as the MATLAB Reachable Set Computation Toolbox CORA. The specific meanings of each parameter can be found in CORA's documentation. please refer to the [manual of CORA](result_picture/Cora2021Manual.pdf).**
```cpp
    NonlinearSys<double> mysys(f, 2, 0, 2);
    ReachOptions<double> options;

    //create R0
    Vector_t<double> center(2);
    Vector_t<double> c(2);
    c << 1.4, 2.4;
    Matrix_t<double> generators(2,1);
    Matrix_t<double> G(2,2);
    G<< 0.17,0,
                 0,0.06;
    Zonotope<double> R0_(c,G);

    center << 1.4, 2.46;
    generators<< 0.17,
                 0;

    options.set_R0(R0_);

    options.set_time_step(0.005);
    options.set_taylor_terms(4);
    options.set_zonotope_order(50);
    options.set_intermediate_order(50);
    options.set_error_order(20);
    options.set_alg("lin");
    options.set_tensor_order(3);

    options.set_tFinal(6.74);
    options.set_tStart(0);

    options.set_usekrylovError(1);
    options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));
```
### 2.2.4 Invoking the Boundary-based Method for Computing Outer-approximations of Reachable Sets
This step invokes our boundary-based method for computing outer-approximations of reachable sets. Please refer to **Section 2.1.1** for the meanings of various parameters.
```cpp
vector<ReachableSet<double>> BdReachset = OverApprox::BdReach(mysys, options, R0_);
```
### 2.2.5 The Plotting of Results
For plotting the graphical results, we utilize the lightweight plotting library **Matplotlib for C++**." For specific usage instructions,please refer to [Matplotlib for C++ Documentation](https://matplotlib-cpp.readthedocs.io/en/latest/index.html).
```cpp
plt::figure_size(1200, 780);
for(int i = 0; i < BdReachset.size(); i++){
    Plotter::plotReach(BdReachset[i], 1, 2, "b");
}
plt::show();
```
### 2.2.6 Results Display
**We employ both the BdryReach and CORA methods to compute the outer-approximation of the reachable set starting from the initial region ([1.23, 1.57], [2.34, 2.46]) over the time interval 0 to 6.74 seconds. The blue region represents the results obtained by the BdryReach method, while the red region corresponds to the results from CORA Computations. It is evident that the outer-approximation computed by BdryReach exhibits significantly higher accuracy compared to CORA.**
<p align="center">
  <img src=result_picture/2.2.6.png>
</p>


## 3 Getting Started

### 3.1 Tool Installation


#### 3.1.1 Eigen3

```bash
sudo apt-get install libeigen3-dev
```

#### 3.1.2 Python

```bash
sudo apt-get install python-dev
```

#### 3.1.3 Capd

```bash
sudo apt install libtool
autoreconf --install
mkdir capd_nogui
cd capd_nogui
../capd/configure --prefix=/usr/local --without-gui --without-mpfr
make
sudo make install
```


### 3.1.4 BdryReach Toolkit Installation and Compilation of Test Cases

```bash
cd BdryReach/
mkdir build
cd build
cmake ..
make
```

## 3.2 Load docker

## 3.3 Simple user Guide 

**In this subsection we will show how to compute an inner-approximation using our tool BdryReach. The main source code is located in ./BdryReach_code/include/underApprox/underApprox.h.**

### 3.3.1 Include Files

```cpp
#include <plotter/matplotlibcpp.h>   // Header for computing reachable set outer-approximation
#include <plotter/plotter.h>          // Header for result visualization
#include <underApprox/underApprox.h>  // Header for includes the interface for computing reachable sets under approximation.
```

### 3.3.2 Definition of Differential Equations

**We use the Capd library to define the form of the differential equations. Refer to the Capd documentation on [differential equation systems](https://capd.sourceforge.net/capdDynSys/docs/html/maps.html). Notably, the computation of our method requires validation of the obtained reachable set inner-approximation. Therefore, an additional definition for a time-inverted differential equation is necessary.**

```cpp
double mu = 1.;

void _f(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = -in[1];
    out[1] = -(0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
}

void _fBack(Node/* t*/, Node in[], int /*dimIn*/, Node out[], int/* dimOut*/, Node params[], int noParams){
    out[0] = in[1];
    out[1] = (0.2 - 0.7*sin(in[0]) - 0.05 * in[1]);
}

int dimIn = 3;

int dimOut = 2;
int noParam = 0;
int MaxDerivativeOrder = 3;

IMap f(_f, dimIn,dimOut,noParam,MaxDerivativeOrder);
IMap fBack(_fBack, dimIn,dimOut,noParam,MaxDerivativeOrder);
```

### 3.3.3 Parameter Configuration for Computing Reachable Sets

**We adopt parameter definitions similar to the MATLAB Reachability Analysis Toolbox CORA. For detailed meanings, refer to CORA's documentation.**
```cpp
NonlinearSys<double> mysys(f, 2, 0, 2);
NonlinearSys<double> mysysBack(fBack, 2, 0, 2);

ReachOptions<double> options;

// create R0
Vector_t<double> center(2);
center << 0, 3;
// center << -7.5, 2.8;
Matrix_t<double> generators(2,2);
generators<< 0.1,0,
                0,0.1;

Zonotope<double> R0_(center,generators);

options.set_time_step(0.2);
options.set_taylor_terms(4);
options.set_zonotope_order(50);
options.set_intermediate_order(50);
options.set_error_order(20);
options.set_alg("lin");
options.set_tensor_order(3);

options.set_tFinal(3);
options.set_tStart(0);

options.set_R0(R0_);

options.set_usekrylovError(1);
options.set_max_error(DBL_MAX*Eigen::MatrixXd::Ones(2,1));
```

### 3.3.4 Invoking the Boundary-based Method for Computing the Inner-approximations of Reachable Sets

**This step invokes our boundary-based method for computing inner-approximations of reachable sets. Please refer to the interface for the meanings of various parameters.**
```cpp
vector<Zonotope<double>> underR = UnderApprox::underReachClp(mysys, mysysBack, options, R0_, 3, 1, 0.01, 0.005, 0.005, 50,50);
```
### Inner-approximation of Reachable Set Computation Interface 
```cpp
template <typename Number>
        static vector<Zonotope<Number>> underReachClp(NonlinearSys<Number> mysys, 			
        NonlinearSys<Number> mysysBack, ReachOptions<Number> options, Zonotope<Number> R0, double 	
        overRtime, int steps, double radius, double over_step, double bound_step, int Zover_order)
```
**Parameters:**
* **mysys:** ordinary differential equation for computing reachable sets.
* **mysysBack:** time-inverted ordinary differential equation for result verification.
* **options:** relevant configuration for outer-approximation of reachable set computation in the program.
* **R0:** initial set.
* **overRtime:** step size for inner-approximation of reachable set computation at each step.
* **steps:** number of iterations for inner-approximation of reachable set computation.
* **radius:** maximum allowed length of generator for facets.
* **over_step:** step size for outer-approximation computation for the entire set at each step in inner-approximation computation.
* **bound_step:** step size for outer-approximation computation for the boundary of the set at each step in inner-approximation computation.
* **Zover_order:** limit on the zonotope order for outer-approximation computation for the entire set at each step in inner-approximation computation.
  
### 3.3.5 The Plotting of Results

For plotting the graphical results, we utilize the lightweight plotting library **Matplotlib for C++**." For specific usage instructions,please refer to [Matplotlib for C++ Documentation](https://matplotlib-cpp.readthedocs.io/en/latest/index.html).
```cpp
plt::figure_size(1200, 780);
for(int i = 1; i < underR.size(); i++){
    Plotter::plotZonotope(underR[i], 1, 2, "g");
}
Plotter::plotZonotope(R0_, 1, 2, "k");
plt::show();
```

### 3.3.6 compile and run 

* **./BdryReach_code/CMakeLists.txt** 
```cpp
cmake_minimum_required(VERSION 3.14)
project(reachSolver)

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)

add_subdirectory(src)
# add_subdirectory(Table3_examples)
# add_subdirectory(Table4_examples)
# add_subdirectory(Table5_examples_shortTime)
add_subdirectory(Table5_examples_longTime)


add  add_subdirectory(Table5_examples_longTime) to the last row.

```

* **./BdryReach_code/Table3_examples/CMakeLists.txt** 
```cpp
cmake_minimum_required(VERSION 3.14)
project(reach_solver_test)

set(CMAKE_CXX_STANDARD 14)

set(TARGET_NAME ElectroOsc)
#ElectroOsc
#Roessler
#LotkaVolterra
#Tank6
#BiologicalSystemI
#BiologicalSystemII
#Tank12
find_package(Eigen3 REQUIRED)
find_package(PythonLibs 2.7 REQUIRED QUIET)

add_executable(${TARGET_NAME} ${TARGET_NAME}.cpp)
target_link_libraries(${TARGET_NAME} PUBLIC
        ${PYTHON_LIBRARIES}
        /usr/local/lib/libcapd.so
        Eigen3::Eigen
        glpk
        reach_solver)
```

* **./BdryReach_code** 
```bash
mkdir build
cd build
cmake ..
make
```

### 3.3.7 Results Display and comments

**zonotopic inner-approximation: the computed inner-approximation**
```markdown
[-6.678246781374445 0.07556403903050581 3.861000979026132e-05 -0.2942877976311855 ;
2.822964194048449 0.007304956271076335 0 0.07559828881651579 ;]
time cost: 22.754546
!
time cost: 22.901619   %the computation time of the whole computation.
```
* **Figure display**
<p align="center">
  <img src=result_picture/3.3.7.jpg>
</p>

**The black region represents the initial set , the green region represents the inner-approximation of the reachable set, while, for comparison, the blue region represents the outer-approximation of the reachable set.**
## 3.4 Replicating Experiments 