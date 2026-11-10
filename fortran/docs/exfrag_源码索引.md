# exfrag 源码索引与接口骨架

对应 2026-09-05 源码快照。此附录按源码提取模块、过程签名、显式 INTENT 和直接 CALL；它不是完整 Fortran 语义分析器，不把函数表达式、泛型运算符或间接调用当作直接 CALL。物理含义与隐式状态见[工程说明书](exfrag_%E5%B7%A5%E7%A8%8B%E4%BB%A3%E7%A0%81%E8%AF%B4%E6%98%8E%E4%B9%A6.md)。没有 INTENT 不代表没有输入输出；特别是 anahyp 的 COMMON 状态必须对照主文阅读。

主入口为 anahyp.for:28–3673 的无名主程序；exfrag.for 仅按顺序 INCLUDE GEMMOD.for 和 anahyp.for。

## 模块目录

| 模块 | 源码起止行 | 过程数 |
|---|---|---:|
| `gem_kind` | 4–50 | 0 |
| `trig` | 54–226 | 14 |
| `sierk_stuff` | 231–1082 | 7 |
| `sortsub` | 1086–1519 | 6 |
| `spectra` | 1523–1672 | 7 |
| `light_stuff` | 1676–2150 | 1 |
| `mass_stuff` | 2154–2577 | 6 |
| `functs` | 2581–3046 | 9 |
| `tl_stuff` | 3052–3977 | 12 |
| `stat_mod` | 3981–4112 | 3 |
| `yrast_mod` | 4118–4758 | 8 |
| `angle` | 4764–5955 | 8 |
| `evap` | 5961–6477 | 4 |
| `evaporation` | 6483–7264 | 7 |
| `nonstatistical_mod` | 7270–7975 | 2 |
| `statistical_mod` | 7980–10131 | 16 |
| `gemsub` | 10137–10786 | 7 |

## gem_kind


## trig

### SIND4

[源码：GEMMOD.for:97](../GEMMOD.for#L97)；过程结束行 102。

```fortran
ELEMENTAL FUNCTION SIND4(theta) RESULT(funct_out)
REAL (kind=4),INTENT(IN) :: theta
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### SIND8

[源码：GEMMOD.for:104](../GEMMOD.for#L104)；过程结束行 110。

```fortran
ELEMENTAL FUNCTION SIND8(theta) RESULT(funct_out)
REAL (kind=8),INTENT(IN) :: theta
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### COSD4

[源码：GEMMOD.for:113](../GEMMOD.for#L113)；过程结束行 119。

```fortran
ELEMENTAL FUNCTION COSD4(theta) RESULT(funct_out)
REAL (kind=4),INTENT(IN) :: theta
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### COSD8

[源码：GEMMOD.for:121](../GEMMOD.for#L121)；过程结束行 127。

```fortran
ELEMENTAL FUNCTION COSD8(theta) RESULT(funct_out)
REAL (kind=8),INTENT(IN) :: theta
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### TAND4

[源码：GEMMOD.for:130](../GEMMOD.for#L130)；过程结束行 137。

```fortran
ELEMENTAL  FUNCTION TAND4(theta) RESULT(funct_out)
REAL (kind=4), INTENT(IN) :: theta
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### TAND8

[源码：GEMMOD.for:139](../GEMMOD.for#L139)；过程结束行 146。

```fortran
ELEMENTAL  FUNCTION TAND8(theta) RESULT(funct_out)
REAL (kind=8), INTENT(IN) :: theta
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### ASIND4

[源码：GEMMOD.for:148](../GEMMOD.for#L148)；过程结束行 155。

```fortran
ELEMENTAL FUNCTION ASIND4(x) RESULT(funct_out)
REAL (kind=4), INTENT(IN) :: x
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### ASIND8

[源码：GEMMOD.for:157](../GEMMOD.for#L157)；过程结束行 164。

```fortran
ELEMENTAL FUNCTION ASIND8(x) RESULT(funct_out)
REAL (kind=8), INTENT(IN) :: x
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### ACOSD4

[源码：GEMMOD.for:167](../GEMMOD.for#L167)；过程结束行 174。

```fortran
ELEMENTAL FUNCTION ACOSD4(x) RESULT(funct_out)
REAL (kind=4), INTENT(IN) :: x
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### ACOSD8

[源码：GEMMOD.for:176](../GEMMOD.for#L176)；过程结束行 183。

```fortran
ELEMENTAL FUNCTION ACOSD8(x) RESULT(funct_out)
REAL (kind=8), INTENT(IN) :: x
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### ATAND4

[源码：GEMMOD.for:185](../GEMMOD.for#L185)；过程结束行 192。

```fortran
ELEMENTAL FUNCTION ATAND4(x) RESULT(funct_out)
REAL (kind=4), INTENT(IN) :: x
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### ATAND8

[源码：GEMMOD.for:194](../GEMMOD.for#L194)；过程结束行 201。

```fortran
ELEMENTAL FUNCTION ATAND8(x) RESULT(funct_out)
REAL (kind=8), INTENT(IN) :: x
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### ATAN2D4

[源码：GEMMOD.for:203](../GEMMOD.for#L203)；过程结束行 213。

```fortran
ELEMENTAL FUNCTION ATAN2D4(x,y) RESULT(funct_out)
REAL (kind=4), INTENT(IN) :: x,y
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### ATAN2D8

[源码：GEMMOD.for:215](../GEMMOD.for#L215)；过程结束行 225。

```fortran
ELEMENTAL FUNCTION ATAN2D8(x,y) RESULT(funct_out)
REAL (kind=8), INTENT(IN) :: x,y
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## sierk_stuff

### jmax_sierk

[源码：GEMMOD.for:638](../GEMMOD.for#L638)；过程结束行 700。

```fortran
SUBROUTINE jmax_sierk (z_in,a_in,Jmax_s)
REAL (KIND=r4), INTENT(IN) :: Z_in, A_in
REAL (KIND=r4), INTENT(OUT), OPTIONAL :: Jmax_s
```

直接 CALL：`lpoly`。

### bfis_sierk

[源码：GEMMOD.for:702](../GEMMOD.for#L702)；过程结束行 761。

```fortran
SUBROUTINE bfis_sierk(j,bfis_s)
REAL (KIND=r4), INTENT(IN) :: J
REAL (KIND=r4), INTENT(OUT) :: bfis_s
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### yrast_sierk

[源码：GEMMOD.for:763](../GEMMOD.for#L763)；过程结束行 796。

```fortran
SUBROUTINE yrast_sierk(J,E_rot_s)
REAL (KIND=r4), INTENT(IN) :: J
REAL (KIND=r4), INTENT(OUT) :: E_rot_s
```

直接 CALL：`lpoly`。

### M_inertia_sierk

[源码：GEMMOD.for:798](../GEMMOD.for#L798)；过程结束行 934。

```fortran
SUBROUTINE M_inertia_sierk(J,M_inertia_1,M_inertia_2,M_inertia_3)
REAL (KIND=r4), INTENT(IN) :: J
REAL (KIND=r4), INTENT(OUT), OPTIONAL :: m_inertia_1,m_inertia_2, m_inertia_3
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### bs_sierk

[源码：GEMMOD.for:936](../GEMMOD.for#L936)；过程结束行 958。

```fortran
SUBROUTINE bs_sierk (j,bs_s)
REAL (kind=r4), INTENT (IN) :: j
REAL (kind=r4), INTENT (OUT) :: bs_s
```

直接 CALL：`lpoly`。

### lpoly

[源码：GEMMOD.for:960](../GEMMOD.for#L960)；过程结束行 991。

```fortran
SUBROUTINE lpoly (x,n,pl)
INTEGER, INTENT(IN) :: n
REAL (KIND=r8), INTENT(IN)  :: x
REAL (KIND=r8), DIMENSION(:), INTENT(OUT) :: pl
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### fisrot

[源码：GEMMOD.for:993](../GEMMOD.for#L993)；过程结束行 1081。

```fortran
SUBROUTINE fisrot (Z,A,al,bf)
REAL (KIND=r4), INTENT(IN) :: a,z,al
REAL (KIND=r4), INTENT(OUT) :: bf
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## sortsub

### write_file_info

[源码：GEMMOD.for:1101](../GEMMOD.for#L1101)；过程结束行 1338。

```fortran
SUBROUTINE write_file_info(unit_out)
INTEGER, INTENT(IN) :: unit_out
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### get_event

[源码：GEMMOD.for:1340](../GEMMOD.for#L1340)；过程结束行 1416。

```fortran
SUBROUTINE get_event (l,I_read_error)
INTEGER, INTENT(OUT) :: l
INTEGER, INTENT(OUT) :: I_read_error
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### read_file_header

[源码：GEMMOD.for:1418](../GEMMOD.for#L1418)；过程结束行 1429。

```fortran
SUBROUTINE read_file_header(unit_in,I_read_error)
INTEGER, INTENT(IN) :: unit_in
INTEGER, INTENT(OUT) :: i_read_error
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### open_event_file

[源码：GEMMOD.for:1431](../GEMMOD.for#L1431)；过程结束行 1474。

```fortran
SUBROUTINE open_event_file(unit_in,I_read_error,unit_out)
INTEGER, INTENT(IN) :: unit_in,unit_out
INTEGER, INTENT(OUT) :: I_read_error
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### start_sort

[源码：GEMMOD.for:1476](../GEMMOD.for#L1476)；过程结束行 1492。

```fortran
SUBROUTINE start_sort(unit_in,unit_out,I_read_error)
INTEGER, INTENT(IN) :: unit_in,unit_out
INTEGER, INTENT(OUT) :: I_read_error
```

直接 CALL：`open_event_file`、`read_file_header`、`write_file_info`。

### name

[源码：GEMMOD.for:1494](../GEMMOD.for#L1494)；过程结束行 1517。

```fortran
FUNCTION name(Z) RESULT(funct_out)
INTEGER, INTENT(IN) :: z
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## spectra

### hist_ek

[源码：GEMMOD.for:1541](../GEMMOD.for#L1541)；过程结束行 1582。

```fortran
SUBROUTINE hist_ek(ek,de,mode,yield)
REAL (KIND=r4), INTENT(IN) :: ek,de,yield
INTEGER, INTENT(IN) :: mode
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### hist_normalize

[源码：GEMMOD.for:1585](../GEMMOD.for#L1585)；过程结束行 1597。

```fortran
SUBROUTINE hist_normalize(gamma_total,j)
REAL (KIND=r4), INTENT(IN) :: gamma_total
INTEGER, INTENT(IN) :: j
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### hist_zero

[源码：GEMMOD.for:1600](../GEMMOD.for#L1600)；过程结束行 1603。

```fortran
SUBROUTINE hist_zero
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### hist_zero_acc

[源码：GEMMOD.for:1606](../GEMMOD.for#L1606)；过程结束行 1609。

```fortran
SUBROUTINE hist_zero_acc
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### hist_write

[源码：GEMMOD.for:1611](../GEMMOD.for#L1611)；过程结束行 1621。

```fortran
SUBROUTINE hist_write(filename)
CHARACTER (LEN=44), INTENT(IN) :: filename
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### hist_read

[源码：GEMMOD.for:1624](../GEMMOD.for#L1624)；过程结束行 1634。

```fortran
SUBROUTINE hist_read(filename)
CHARACTER (LEN=44), INTENT(IN) :: filename
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### hist_print

[源码：GEMMOD.for:1636](../GEMMOD.for#L1636)；过程结束行 1671。

```fortran
SUBROUTINE hist_print(plb,lmin,lmax,dif)
real(kind=r4), INTENT(IN) :: plb
REAL (KIND=r4), INTENT(IN):: lmax
REAL (KIND=r4), INTENT(IN) :: dif
INTEGER, INTENT(IN) :: lmin
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## light_stuff

### write_mult

[源码：GEMMOD.for:2051](../GEMMOD.for#L2051)；过程结束行 2149。

```fortran
SUBROUTINE write_mult (num_residue)
REAL (KIND=r4), INTENT(IN) :: num_residue
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## mass_stuff

### getmass

[源码：GEMMOD.for:2187](../GEMMOD.for#L2187)；过程结束行 2244。

```fortran
SUBROUTINE getmass(iZ,iA,mass_exp,mass_ld,beta_two)
INTEGER, INTENT(IN) :: iZ,iA
REAL (KIND=r4), INTENT(OUT), OPTIONAL :: mass_ld,mass_exp,beta_two
```

直接 CALL：`massfr`。

### prepare_mass

[源码：GEMMOD.for:2246](../GEMMOD.for#L2246)；过程结束行 2296。

```fortran
SUBROUTINE prepare_mass
```

直接 CALL：`prepare_mass_tf`、`get_pairing`。

### prepare_mass_tf

[源码：GEMMOD.for:2298](../GEMMOD.for#L2298)；过程结束行 2374。

```fortran
SUBROUTINE prepare_mass_tf
```

直接 CALL：`get_pairing`。

### massfr

[源码：GEMMOD.for:2376](../GEMMOD.for#L2376)；过程结束行 2478。

```fortran
PURE SUBROUTINE massfr(Z,AA,mass)
REAL (KIND=r4), INTENT(IN) :: Z,AA
REAL (KIND=r4), INTENT(OUT) :: mass
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### droplet_mass

[源码：GEMMOD.for:2480](../GEMMOD.for#L2480)；过程结束行 2537。

```fortran
PURE SUBROUTINE droplet_mass(Z,A,mass)
REAL (KIND=r4), INTENT(IN) :: Z,A
REAL (KIND=r4), INTENT(OUT) :: mass
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### get_pairing

[源码：GEMMOD.for:2539](../GEMMOD.for#L2539)；过程结束行 2575。

```fortran
PURE SUBROUTINE get_pairing(iz,ia,pairing)
INTEGER, INTENT(IN) :: iz,ia
REAL (KIND=r4), INTENT(OUT) :: pairing
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## functs

### gauss

[源码：GEMMOD.for:2599](../GEMMOD.for#L2599)；过程结束行 2626。

```fortran
FUNCTION gauss () RESULT (FUNCT_out)
```

直接 CALL：`random_number`。

### cl_go

[源码：GEMMOD.for:2628](../GEMMOD.for#L2628)；过程结束行 2687。

```fortran
PURE FUNCTION cl_go(J1,J2,J,M1,M2,M) RESULT(funct_out)
INTEGER, INTENT(IN) :: J1,J2,J,M1,M2,M
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### factorial

[源码：GEMMOD.for:2689](../GEMMOD.for#L2689)；过程结束行 2698。

```fortran
SUBROUTINE factorial ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### exp_decay

[源码：GEMMOD.for:2700](../GEMMOD.for#L2700)；过程结束行 2710。

```fortran
FUNCTION exp_decay(decay_width) RESULT(funct_out)
REAL (kind=r4), INTENT(IN) :: decay_width
```

直接 CALL：`random_number`。

### transform

[源码：GEMMOD.for:2712](../GEMMOD.for#L2712)；过程结束行 2777。

```fortran
PURE SUBROUTINE transform(theta1,phi1,theta2,phi2,theta3,phi3)
REAL (KIND=r4), INTENT(IN) :: theta2,phi2
REAL (KIND=r4), INTENT(IN) :: theta1,phi1
REAL (KIND=r4), INTENT(OUT) :: theta3,phi3
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### addvectors

[源码：GEMMOD.for:2779](../GEMMOD.for#L2779)；过程结束行 2820。

```fortran
PURE FUNCTION addvectors (v1,v2) RESULT(v3)
TYPE (vell), INTENT(IN) :: v1,v2
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### ran_j

[源码：GEMMOD.for:2822](../GEMMOD.for#L2822)；过程结束行 2871。

```fortran
SUBROUTINE ran_j(x,ac,U,aden,rua,xmax)
REAL (KIND=r4), INTENT(OUT) :: x
REAL (KIND=r4), INTENT(IN) :: aden,rua,xmax,U,ac
```

直接 CALL：`random_number`。

### get_little_a

[源码：GEMMOD.for:2873](../GEMMOD.for#L2873)；过程结束行 3037。

```fortran
SUBROUTINE get_little_a(A,Z,J,U,surf,shell_corr,E_rot,aden, entropy, pair_energy,aa_scale)
REAL (KIND=r4), INTENT(IN) :: A,Z,J,surf,aa_scale,pair_energy, shell_corr,E_rot
REAL (KIND=r4), INTENT(INOUT) :: U
REAL (KIND=r4), INTENT(OUT) :: aden,entropy
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### zbeta

[源码：GEMMOD.for:3039](../GEMMOD.for#L3039)；过程结束行 3044。

```fortran
ELEMENTAL FUNCTION zbeta(A) RESULT(funct_out)
REAL (KIND=r4), INTENT(IN) :: A
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## tl_stuff

### get_tl

[源码：GEMMOD.for:3124](../GEMMOD.for#L3124)；过程结束行 3296。

```fortran
SUBROUTINE get_tl(l,ek,tl,tl_up,tl_down)
INTEGER, INTENT(IN) :: l
REAL (KIND=r4), INTENT(IN) :: ek
REAL (KIND=r4), DIMENSION(0:l_max_quantum), INTENT(OUT) :: tl
REAL (KIND=r4), INTENT(OUT) :: tl_up,tl_down
```

直接 CALL：`get_tll`。

### get_tl_half

[源码：GEMMOD.for:3298](../GEMMOD.for#L3298)；过程结束行 3335。

```fortran
SUBROUTINE get_tl_half(mode,iz2,l_min,ek_min)
INTEGER, INTENT(IN) :: mode,iz2, l_min
REAL   (KIND=r4), INTENT(OUT) :: ek_min
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### prepare_deformed_emission

[源码：GEMMOD.for:3338](../GEMMOD.for#L3338)；过程结束行 3425。

```fortran
SUBROUTINE prepare_deformed_emission (ratio,expand)
REAL (KIND=r4), INTENT(IN) :: ratio,expand
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### read_tl

[源码：GEMMOD.for:3427](../GEMMOD.for#L3427)；过程结束行 3563。

```fortran
SUBROUTINE read_tl ()
```

直接 CALL：`prepare_interpolation`、`prepare_plm2`。

### clear_look_up

[源码：GEMMOD.for:3565](../GEMMOD.for#L3565)；过程结束行 3570。

```fortran
SUBROUTINE clear_look_up ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### prepare_plm2

[源码：GEMMOD.for:3572](../GEMMOD.for#L3572)；过程结束行 3653。

```fortran
SUBROUTINE prepare_plm2 ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### prepare_interpolation

[源码：GEMMOD.for:3655](../GEMMOD.for#L3655)；过程结束行 3754。

```fortran
SUBROUTINE prepare_interpolation ()
```

直接 CALL：`inverse`。

### get_tll

[源码：GEMMOD.for:3756](../GEMMOD.for#L3756)；过程结束行 3795。

```fortran
PURE SUBROUTINE get_tll(c,ek,eetl)
REAL (KIND=r4), DIMENSION(3), INTENT(IN) :: c
REAL (KIND=r4), INTENT(IN) :: ek
REAL (KIND=r4), INTENT(OUT) :: eetl
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### inverse

[源码：GEMMOD.for:3797](../GEMMOD.for#L3797)；过程结束行 3818。

```fortran
PURE SUBROUTINE inverse(m,m_inv)
REAL (KIND=r4), INTENT(IN), DIMENSION(:,:)  :: m
REAL (KIND=r4), INTENT(OUT), DIMENSION(:,:) :: m_inv
```

直接 CALL：`sle`。

### sle

[源码：GEMMOD.for:3821](../GEMMOD.for#L3821)；过程结束行 3844。

```fortran
PURE SUBROUTINE sle(m,y)
REAL (KIND=r4), INTENT(INOUT), DIMENSION(:,:) :: m
REAL (KIND=r4), INTENT(INOUT), DIMENSION(:)   :: y
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### prepare_tl

[源码：GEMMOD.for:3846](../GEMMOD.for#L3846)；过程结束行 3947。

```fortran
SUBROUTINE prepare_tl(mode,iz2,ia2)
INTEGER, INTENT(IN) :: mode,iz2,ia2
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### DE_ARRAY_tl

[源码：GEMMOD.for:3950](../GEMMOD.for#L3950)；过程结束行 3975。

```fortran
SUBROUTINE DE_ARRAY_tl
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## stat_mod

### alloc_storage

[源码：GEMMOD.for:4025](../GEMMOD.for#L4025)；过程结束行 4032。

```fortran
SUBROUTINE alloc_storage ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### get_level_density

[源码：GEMMOD.for:4034](../GEMMOD.for#L4034)；过程结束行 4078。

```fortran
SUBROUTINE get_level_density(U,J,aden,M_Inertia, entropy,level_density)
REAL (KIND=r4), INTENT(IN) :: U,J
REAL (KIND=r4), INTENT(IN) :: aden,M_inertia,entropy
REAL (KIND=r4), INTENT(OUT) :: level_density
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### collective_enhancement

[源码：GEMMOD.for:4081](../GEMMOD.for#L4081)；过程结束行 4110。

```fortran
SUBROUTINE collective_enhancement(IZ,IA,U,aden,inertia, level_density)
INTEGER, INTENT(IN) :: iz,ia
REAL (KIND=r4), INTENT(IN) :: inertia,U,aden
REAL (KIND=r4), INTENT(INOUT) :: level_density
```

直接 CALL：`getmass`。


## yrast_mod

### sadfit

[源码：GEMMOD.for:4136](../GEMMOD.for#L4136)；过程结束行 4357。

```fortran
SUBROUTINE sadfit(Z0,A0,J0)
REAL (KIND=r4), INTENT(IN) :: Z0,A0,J0
```

直接 CALL：`jmax_sierk`、`bfis_sierk`、`getmass`。

### read_sad

[源码：GEMMOD.for:4360](../GEMMOD.for#L4360)；过程结束行 4396。

```fortran
SUBROUTINE read_sad ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### de_c

[源码：GEMMOD.for:4399](../GEMMOD.for#L4399)；过程结束行 4402。

```fortran
subroutine de_c
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### yrastld

[源码：GEMMOD.for:4405](../GEMMOD.for#L4405)；过程结束行 4542。

```fortran
SUBROUTINE yrastld (Z,A,l,E_rot)
REAL (KIND=r4), INTENT(IN) :: Z,A,l
REAL (KIND=r4), INTENT(OUT) :: E_rot
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### yrast

[源码：GEMMOD.for:4544](../GEMMOD.for#L4544)；过程结束行 4639。

```fortran
SUBROUTINE yrast(Z,A,J,K,ratio,E_rot)
REAL (KIND=r4), INTENT(IN) :: Z,A,J,K,ratio
REAL (KIND=r4), INTENT(OUT) :: E_rot
```

直接 CALL：`yrastld`、`jmax_sierk`、`yrast_sierk`、`m_inertia_sierk`。

### cubic

[源码：GEMMOD.for:4641](../GEMMOD.for#L4641)；过程结束行 4659。

```fortran
PURE FUNCTION cubic(a,b,c,d,e,f) RESULT(funct_out)
REAL (KIND=r4), INTENT(IN) :: a,b,c,d,e,f
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### saddle

[源码：GEMMOD.for:4662](../GEMMOD.for#L4662)；过程结束行 4727。

```fortran
SUBROUTINE saddle(Z1,A1,Z2,A2,E_bar,m_inertia_eff)
REAL (KIND=r4), INTENT(IN) :: A1,A2,Z1,Z2
REAL (KIND=r4), INTENT(OUT) :: E_bar,m_inertia_eff
```

直接 CALL：`getmass`。

### def_rot

[源码：GEMMOD.for:4729](../GEMMOD.for#L4729)；过程结束行 4756。

```fortran
SUBROUTINE def_rot (ratio,bss)
REAL (KIND=r4), INTENT(IN) :: ratio
REAL (KIND=r4), INTENT(OUT) :: bss
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## angle

### sphere_imf

[源码：GEMMOD.for:4784](../GEMMOD.for#L4784)；过程结束行 5157。

```fortran
SUBROUTINE sphere_imf (parent,daughter,U_thermal)
REAL, INTENT(IN) :: u_thermal
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
```

直接 CALL：`get_little_a`、`ran_j`、`random_number`、`transform`。

### sphere_lp

[源码：GEMMOD.for:5160](../GEMMOD.for#L5160)；过程结束行 5323。

```fortran
SUBROUTINE sphere_lp(parent,daughter,ll,l_p_s,e_kinetic)
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
INTEGER, INTENT(IN) :: ll
REAL (KIND=r4), INTENT(IN) :: l_P_s,e_kinetic
```

直接 CALL：`random_number`、`transform`、`find_theta`。

### quantum_lp

[源码：GEMMOD.for:5326](../GEMMOD.for#L5326)；过程结束行 5446。

```fortran
SUBROUTINE quantum_lp(parent,daughter,ll,l_p_s,e_kinetic)
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
INTEGER, INTENT(IN) :: ll
```

直接 CALL：`random_number`、`find_theta`。

### quantum_ang_generator

[源码：GEMMOD.for:5449](../GEMMOD.for#L5449)；过程结束行 5553。

```fortran
SUBROUTINE quantum_ang_generator ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### find_theta

[源码：GEMMOD.for:5555](../GEMMOD.for#L5555)；过程结束行 5581。

```fortran
SUBROUTINE find_theta(l,m,angle)
INTEGER, INTENT(IN) ::l,m
REAL (KIND=r4), INTENT(OUT) :: angle
```

直接 CALL：`random_number`。

### quantum_imf

[源码：GEMMOD.for:5584](../GEMMOD.for#L5584)；过程结束行 5861。

```fortran
SUBROUTINE quantum_imf (parent,daughter,U_saddle,e_kinetic, U_scission)
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
REAL (KIND=r4) , INTENT(IN) :: U_saddle,e_kinetic
REAL (KIND=r4) , INTENT(IN) :: U_scission
```

直接 CALL：`m_inertia_sierk`、`bs_sierk`、`get_little_a`、`random_number`、`jmax_sierk`、`yrast_sierk`。

### fast_fission_angle

[源码：GEMMOD.for:5863](../GEMMOD.for#L5863)；过程结束行 5940。

```fortran
SUBROUTINE fast_fission_angle(parent,daughter,ek,U_scission)
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
REAL (KIND=r4), INTENT(IN) :: ek,u_scission
```

直接 CALL：`random_number`、`jmax_sierk`、`yrast_sierk`、`get_little_a`。

### DE_ARRAY_angle

[源码：GEMMOD.for:5943](../GEMMOD.for#L5943)；过程结束行 5952。

```fortran
SUBROUTINE DE_ARRAY_angle
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## evap

### evaporation_simple

[源码：GEMMOD.for:6063](../GEMMOD.for#L6063)；过程结束行 6171。

```fortran
SUBROUTINE evaporation_simple(mode,Z0,A0,Ex0,J0)
INTEGER, INTENT(IN) :: mode
REAL (KIND=r4), INTENT(IN) :: Z0,A0,EX0,J0
```

直接 CALL：`get_pairing`、`getmass`、`yrast`、`get_little_a`、`light_particle_quick`、`light_particle`。

### light_particle

[源码：GEMMOD.for:6173](../GEMMOD.for#L6173)；过程结束行 6271。

```fortran
SUBROUTINE light_particle (mode,Ex0)
INTEGER, INTENT(IN) :: mode
REAL (KIND=r4), INTENT(IN) :: EX0
```

直接 CALL：`yrast`、`get_little_a`、`get_level_density`、`random_number`。

### light_particle_quick

[源码：GEMMOD.for:6273](../GEMMOD.for#L6273)；过程结束行 6460。

```fortran
SUBROUTINE light_particle_quick (mode,Ex0)
INTEGER, INTENT(IN) :: mode
REAL (KIND=r4), INTENT(IN) :: EX0
```

直接 CALL：`yrast`、`get_little_a`、`get_level_density`、`random_number`。

### prepare_evap

[源码：GEMMOD.for:6463](../GEMMOD.for#L6463)；过程结束行 6476。

```fortran
SUBROUTINE prepare_evap ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## evaporation

### Hauser_Feshbach

[源码：GEMMOD.for:6532](../GEMMOD.for#L6532)；过程结束行 6832。

```fortran
SUBROUTINE Hauser_Feshbach(mode,z0,a0,Ex0,J0,k0,ratio)
INTEGER, INTENT(IN) :: mode
REAL (KIND=r4), INTENT(IN)  :: J0,K0,Z0,A0,Ex0,ratio
```

直接 CALL：`get_pairing`、`getmass`、`prepare_tl`、`clear_look_up`、`random_number`、`s2_loop`、`yrast`、`ell_loop`。

### s2_loop

[源码：GEMMOD.for:6835](../GEMMOD.for#L6835)；过程结束行 6895。

```fortran
SUBROUTINE s2_loop (mode,z0,a0,Ex0,k0,ratio)
INTEGER, INTENT(IN) :: mode
REAL (KIND=r4), INTENT(IN)  :: K0,Z0,A0,Ex0,ratio
```

直接 CALL：`k2_loop`、`yrast`、`ek_loop`。

### ek_loop

[源码：GEMMOD.for:6897](../GEMMOD.for#L6897)；过程结束行 6985。

```fortran
SUBROUTINE ek_loop (mode,z0,a0,Ex0,k0)
INTEGER, INTENT(IN) :: mode
REAL (KIND=r4), INTENT(IN)  :: K0,Z0,A0,Ex0
```

直接 CALL：`ell_loop`、`get_little_a`、`get_level_density`、`collective_enhancement`。

### ell_loop

[源码：GEMMOD.for:6987](../GEMMOD.for#L6987)；过程结束行 7146。

```fortran
SUBROUTINE ell_loop (mode,k0,j_sub)
INTEGER, INTENT(IN) :: mode
REAL (KIND=r4), INTENT(IN)  :: K0
INTEGER, INTENT(INOUT) :: j_sub
```

直接 CALL：`get_tl`。

### k2_loop

[源码：GEMMOD.for:7148](../GEMMOD.for#L7148)；过程结束行 7223。

```fortran
SUBROUTINE k2_loop (mode,z0,a0,Ex0,k0,ratio)
INTEGER, INTENT(IN) :: mode
REAL (KIND=r4), INTENT(IN)  :: K0,Z0,A0,Ex0,ratio
```

直接 CALL：`yrast`、`ek_loop`。

### prepare_Hauser_Feshbach

[源码：GEMMOD.for:7227](../GEMMOD.for#L7227)；过程结束行 7241。

```fortran
SUBROUTINE prepare_Hauser_Feshbach ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### prepare_5li_5he

[源码：GEMMOD.for:7243](../GEMMOD.for#L7243)；过程结束行 7262。

```fortran
SUBROUTINE prepare_5li_5he ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## nonstatistical_mod

### prepare_nonstatistical

[源码：GEMMOD.for:7407](../GEMMOD.for#L7407)；过程结束行 7428。

```fortran
SUBROUTINE prepare_nonstatistical ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### non_statistical

[源码：GEMMOD.for:7431](../GEMMOD.for#L7431)；过程结束行 7974。

```fortran
SUBROUTINE non_statistical(parent,daughter,n_extra)
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), DIMENSION(3), INTENT(out) :: daughter
INTEGER, INTENT(out) :: n_extra
```

直接 CALL：`prepare_nonstatistical`、`random_number`、`quantum_lp`。


## statistical_mod

### decay

[源码：GEMMOD.for:8015](../GEMMOD.for#L8015)；过程结束行 8112。

```fortran
SUBROUTINE decay(parent,daughter,n_extra,ratio,weight,n_binary)
TYPE (identity), INTENT(INOUT) :: parent
TYPE (identity), INTENT(OUT), DIMENSION(3) :: daughter
INTEGER, INTENT(OUT) :: n_extra
INTEGER, INTENT(INOUT) :: n_binary
REAL (KIND=r4), INTENT(INOUT) :: ratio
REAL (KIND=r4), INTENT(INOUT) :: weight
```

直接 CALL：`non_statistical`、`statistical`、`fission_width`、`viola`、`fast_fission_angle`、`sphere_imf`、`xray`。

### statistical

[源码：GEMMOD.for:8115](../GEMMOD.for#L8115)；过程结束行 8192。

```fortran
SUBROUTINE statistical(parent,daughter,n_extra,ratio,weight)
TYPE (identity), INTENT(INOUT) :: parent
TYPE (identity), INTENT(OUT), DIMENSION(2) :: daughter
INTEGER, INTENT(OUT) :: n_extra
REAL (KIND=r4), INTENT(IN) :: ratio
REAL (KIND=r4), INTENT(INOUT) :: weight
```

直接 CALL：`decay_widths`、`decay_time_wise`、`decay_prob_wise`。

### INIT_decay

[源码：GEMMOD.for:8194](../GEMMOD.for#L8194)；过程结束行 8205。

```fortran
SUBROUTINE INIT_decay ()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### decay_time_wise

[源码：GEMMOD.for:8207](../GEMMOD.for#L8207)；过程结束行 8286。

```fortran
SUBROUTINE decay_time_wise(parent,daughter)
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(OUT),DIMENSION(2) :: daughter
```

直接 CALL：`decay_lp`、`decay_imf`。

### decay_prob_wise

[源码：GEMMOD.for:8288](../GEMMOD.for#L8288)；过程结束行 8510。

```fortran
SUBROUTINE decay_prob_wise(parent,daughter,weight)
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(OUT),DIMENSION(2) :: daughter
```

直接 CALL：`random_number`、`decay_lp`、`decay_imf`。

### decay_imf

[源码：GEMMOD.for:8513](../GEMMOD.for#L8513)；过程结束行 8580。

```fortran
SUBROUTINE decay_imf (parent,daughter,i,decay_time)
INTEGER, INTENT(IN) :: i
REAL (KIND=r4), INTENT(IN) :: decay_time
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
```

直接 CALL：`quantum_imf`、`sphere_imf`。

### decay_lp

[源码：GEMMOD.for:8582](../GEMMOD.for#L8582)；过程结束行 8633。

```fortran
SUBROUTINE decay_lp (parent,daughter,i,decay_time)
INTEGER, INTENT(IN) :: i
REAL (KIND=r4), INTENT(IN) :: decay_time
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
```

直接 CALL：`quantum_lp`、`sphere_lp`。

### decay_widths

[源码：GEMMOD.for:8635](../GEMMOD.for#L8635)；过程结束行 8917。

```fortran
SUBROUTINE decay_widths(parent,rratio,gamma_total)
TYPE (identity), INTENT(INOUT) :: parent
REAL (KIND=r4), INTENT(IN) :: rratio
REAL (KIND=r4), INTENT(out) :: gamma_total
```

直接 CALL：`get_pairing`、`getmass`、`fluctuate`、`def_rot`、`yrast`、`prepare_deformed_emission`、`get_little_a`、`get_level_density`、`collective_enhancement`、`hauser_feshbach`、`evaporation_simple`、`gamma_ray`、`fission_width`、`sadfit`、`imf`、`normalise_fission`。

### gamma_ray

[源码：GEMMOD.for:8919](../GEMMOD.for#L8919)；过程结束行 9047。

```fortran
SUBROUTINE gamma_ray(mode,Z0,A0,Ex0,J0,ratio)
INTEGER, INTENT(IN) :: mode
REAL (KIND=r4), INTENT(IN) :: Z0,A0,Ex0,J0,ratio
```

直接 CALL：`yrast`、`get_little_a`、`get_level_density`、`random_number`。

### IMF

[源码：GEMMOD.for:9050](../GEMMOD.for#L9050)；过程结束行 9120。

```fortran
SUBROUTINE IMF(z0,a0,Ex0,J0)
REAL (KIND=r4), INTENT(IN) :: Z0, A0, Ex0, J0
```

直接 CALL：`complex_frag`。

### Normalise_fission

[源码：GEMMOD.for:9127](../GEMMOD.for#L9127)；过程结束行 9163。

```fortran
SUBROUTINE Normalise_fission(Z0,A0,Ex0,J0)
REAL (KIND=r4), INTENT(IN) :: Z0,A0,Ex0,J0
```

直接 CALL：`fission_width`。

### complex_frag

[源码：GEMMOD.for:9165](../GEMMOD.for#L9165)；过程结束行 9318。

```fortran
SUBROUTINE complex_frag(A0,Z1,A1,Z2,A2,Ex0,J0,kl)
INTEGER, INTENT(OUT) :: kl
REAL (KIND=r4), INTENT(IN) :: A0,Z1,A1,Z2,A2,J0,Ex0
```

直接 CALL：`saddle`、`get_little_a`。

### fission_width

[源码：GEMMOD.for:9320](../GEMMOD.for#L9320)；过程结束行 9555。

```fortran
SUBROUTINE fission_width(Z0,A0,Ex0,J0,fast_fission)
REAL (KIND=r4), INTENT(IN) :: Z0,A0,Ex0,J0
LOGICAL, INTENT(OUT), OPTIONAL :: fast_fission
```

直接 CALL：`yrast`、`jmax_sierk`、`bfis_sierk`、`bs_sierk`、`get_little_a`、`m_inertia_sierk`、`random_number`、`viola`。

### fluctuate

[源码：GEMMOD.for:9558](../GEMMOD.for#L9558)；过程结束行 10051。

```fortran
SUBROUTINE fluctuate(Z0,A0,Ex0,J0,k0,ratio)
REAL (KIND=r4), INTENT(IN) :: A0,Z0,Ex0,j0
REAL (KIND=r4), INTENT(OUT) :: k0,ratio
```

直接 CALL：`yrast`、`jmax_sierk`、`bfis_sierk`、`m_inertia_sierk`、`get_little_a`、`get_level_density`、`random_number`。

### viola

[源码：GEMMOD.for:10054](../GEMMOD.for#L10054)；过程结束行 10072。

```fortran
SUBROUTINE viola(Z,A,E_scission,ek_viola,z1,z2,a1,a2)
INTEGER, INTENT(IN) :: Z,A
REAL (KIND=r4), INTENT(OUT) :: E_scission,ek_viola
INTEGER, INTENT(OUT) :: z1,z2,a1,a2
```

直接 CALL：`getmass`。

### xray

[源码：GEMMOD.for:10074](../GEMMOD.for#L10074)；过程结束行 10129。

```fortran
SUBROUTINE xray(parent,daughter,n_extra)
TYPE (identity), INTENT(IN) :: parent
TYPE (identity), INTENT(INOUT), DIMENSION(3) :: daughter
INTEGER, INTENT(INOUT) :: n_extra
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## gemsub

### INPUT_ANA

[源码：GEMMOD.for:10182](../GEMMOD.for#L10182)；过程结束行 10211。

```fortran
SUBROUTINE INPUT_ANA()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### event

[源码：GEMMOD.for:10214](../GEMMOD.for#L10214)；过程结束行 10284。

```fortran
SUBROUTINE event(cn,weight)
TYPE (identity), INTENT(INOUT) :: cn
REAL (KIND=r4), INTENT(out) :: weight
```

直接 CALL：`init_decay`、`decay_event`。

### decay_event

[源码：GEMMOD.for:10286](../GEMMOD.for#L10286)；过程结束行 10314。

```fortran
RECURSIVE SUBROUTINE decay_event(parent,ratio,weight)
TYPE (identity), INTENT(INOUT) :: parent
REAL (KIND=r4), INTENT(INOUT) :: ratio,weight
```

直接 CALL：`decay`、`decay_event`。

### dmp

[源码：GEMMOD.for:10316](../GEMMOD.for#L10316)；过程结束行 10464。

```fortran
SUBROUTINE dmp(L,weight)
INTEGER, INTENT(IN) :: L
REAL (KIND=r4), INTENT(IN) :: weight
```

直接 CALL：`random_number`。

### fusion

[源码：GEMMOD.for:10466](../GEMMOD.for#L10466)；过程结束行 10573。

```fortran
SUBROUTINE fusion(L,CN)
INTEGER, INTENT(IN) :: L
TYPE (identity), INTENT(out) :: cn
```

直接 CALL：`random_number`。

### input

[源码：GEMMOD.for:10575](../GEMMOD.for#L10575)；过程结束行 10760。

```fortran
SUBROUTINE input(ICZ,ICN,NEXC, ratio_ana)
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### crash

[源码：GEMMOD.for:10762](../GEMMOD.for#L10762)；过程结束行 10784。

```fortran
SUBROUTINE crash()
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。


## external

### ISOTAB

[源码：anahyp.for:3680](../anahyp.for#L3680)；过程结束行 3924。

```fortran
FUNCTION ISOTAB(IZ,IN)
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### sl

[源码：anahyp.for:3926](../anahyp.for#L3926)；过程结束行 3938。

```fortran
function sl(i)
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### xexp

[源码：anahyp.for:3941](../anahyp.for#L3941)；过程结束行 3952。

```fortran
function xexp(x)
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### scou

[源码：anahyp.for:3954](../anahyp.for#L3954)；过程结束行 3968。

```fortran
function scou(ncc, al)
```

直接 CALL：`erf`。

### erf

[源码：anahyp.for:3973](../anahyp.for#L3973)；过程结束行 3986。

```fortran
subroutine erf(x,er)
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

### LORENTZ

[源码：anahyp.for:3991](../anahyp.for#L3991)；过程结束行 4009。

```fortran
subroutine LORENTZ(BX, BY, BZ, A0, AX, AY, AZ, AP0, APX, APY, APZ)
```

无直接 CALL（函数表达式及全局状态访问仍可能存在）。

