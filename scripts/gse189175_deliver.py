# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import numpy as np,json,os
os.environ['MPLCONFIGDIR']=str(Path('work/mplcache').resolve())
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=json.loads(Path('work/GSE189175_score_results.json').read_text());audit=json.loads(Path('work/GSE189175_cluster_audit.json').read_text())
z=np.load('work/GSE189175_cluster_inputs.npz');s=np.load('work/GSE189175_scores.npz')
core=r['annotation']['core_clusters'];mixed=r['annotation']['mixed_sensitivity_clusters']
marker_groups={k:[g for g in v if g in audit['clusters'][0]['markers']] for k,v in audit['marker_sets'].items()};gs=sum(marker_groups.values(),[])
fig,ax=plt.subplots(figsize=(15,8))
for yi,row in enumerate(audit['clusters']):
 for xi,g in enumerate(gs):
  v=row['markers'][g];maximum=max(a['markers'][g]['mean'] for a in audit['clusters'])
  color=v['mean']/maximum if maximum else 0
  sc=ax.scatter(xi,yi,s=160*v['detection'],c=[color],cmap='viridis',vmin=0,vmax=1,edgecolors='none')
ax.set_xticks(range(len(gs)),gs,rotation=70,ha='right');ax.set_yticks(range(18),[f"C{i} (n={audit['clusters'][i]['n']:,})" for i in range(18)])
for i,label in enumerate(ax.get_yticklabels()):
 if i in core:label.set_color('#12683b');label.set_fontweight('bold')
 elif i in mixed:label.set_color('#a36200')
ax.invert_yaxis();ax.set_ylim(17.6,-.7);ax.set_xlim(-.7,len(gs)-.3)
start=0
for group,markers in marker_groups.items():
 if start:ax.axvline(start-.5,color='#bbbbbb',lw=.7)
 start+=len(markers)
ax.set_title('GSE189175: lineage marker audit before repair-score contrasts',loc='left',pad=15)
fig.colorbar(sc,ax=ax,shrink=.65,pad=.02,label='Mean log expression / maximum across clusters for that gene')
ax.grid(axis='y',color='#eeeeee');ax.set_axisbelow(True)
fig.text(.05,.02,'Dot area: fraction of nuclei with detected transcript. Green: core hepatic clusters. Amber: mixed lineage sensitivity. Labels inferred; malignancy not assigned.',fontsize=9)
fig.tight_layout(rect=(0,.07,1,1))
fig.savefig('outputs/GSE189175_marker_audit.png',dpi=180);fig.savefig('outputs/GSE189175_marker_audit.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4.5))
for i,row in enumerate(r['patients']):
 a=row['selections']['core'];vals=[a['raw'],a['models']['cycle']['difference'],a['models']['identity']['difference'],a['models']['zone']['difference']]
 ax.plot(range(4),vals,marker='o',lw=1.7,label='Patient '+row['patient'])
ax.axhline(0,color='black',lw=.8,linestyle='--');ax.set_xticks(range(4),['Unadjusted','Depth + S/G2M','+ Hepatic identity','+ Zonation contrast'])
ax.set_ylabel('Repair-score difference (tumor minus adjacent)');ax.set_title('GSE189175: descriptive sensitivity across three patients',loc='left')
ax.legend(frameon=False);ax.grid(axis='y',alpha=.2);fig.text(.08,.025,'No confidence intervals or significance tests. Reannotated hepatic nuclei; n = 3 independent patients.',fontsize=8)
fig.tight_layout(rect=(0,.07,1,1));fig.savefig('outputs/GSE189175_control_sensitivity.png',dpi=200);fig.savefig('outputs/GSE189175_control_sensitivity.pdf');plt.close(fig)
# Independent R QR check of weighted fits.
C=np.column_stack([np.log1p(z['total']),s['S'],s['G2M'],s['identity'],s['portal']-s['central']])
mask=s['coremask'];t=(z['tissue']=='Tumor').astype(int)
np.savetxt('work/GSE189175_qr_check.txt',np.column_stack([z['patients'][mask].astype(int),t[mask],s['repair'][mask],C[mask]]),fmt='%.12g')
md=['# GSE189175: bağımsız kohortta betimleyici sınama','','10 Ekim 2026','','## Sonuç','','Yeniden anotasyonla seçilen 25.977 hepatosit benzeri çekirdekte üç hasta çifti incelendi. Ham farklar 2 negatif, 1 pozitif; genişletilmiş kontrol modeli de 2 negatif, 1 pozitif. Ancak pozitif hastanın sonucu model ve anotasyon seçimine duyarlı, kovaryat örtüşmesi çok sınırlı. İlk kohorttaki tüm hastalarda ortak artış bulunmamasıyla uyumlu bir heterojenlik görüldü; bu, doğrulanmış mekanizma veya kesin dış replikasyon iddiası değildir.','','## Ana hasta sonuçları','','Pozitif değer tümör dokusunda daha yüksek 47 genlik onarım imzası puanını gösterir. Ölçek ortalama log-normalize ifadedir; yüzde değişim değildir.','','| Hasta | Komşu doku çekirdekleri | Tümör çekirdekleri | Ham fark | Derinlik + S/G2M | + Kimlik | + Zonasyon |','|---|---:|---:|---:|---:|---:|---:|']
for row in r['patients']:
 a=row['selections']['core'];md.append(f"| {row['patient']} | {a['n_N']} | {a['n_T']} | {a['raw']:+.4f} | {a['models']['cycle']['difference']:+.4f} | {a['models']['identity']['difference']:+.4f} | {a['models']['zone']['difference']:+.4f} |")
md+=['','## Anotasyon ve örtüşme duyarlılığı','','| Hasta | Ana kontrollü fark | Karışık kümeler dahil | Katı belirteç alt kümesi | Ortak aralık N / T | Ortak aralık kontrollü fark |','|---|---:|---:|---:|---:|---:|']
for row in r['patients']:
 a=row['selections']['core'];b=row['selections']['mixed_included'];c=row['selections']['strict_markers']
 common=f"{a['common_range_models']['zone']['difference']:+.4f}" if 'common_range_models' in a else 'hesaplanmadı'
 md.append(f"| {row['patient']} | {a['models']['zone']['difference']:+.4f} | {b['models']['zone']['difference']:+.4f} | {c['models']['zone']['difference']:+.4f} | {a['common_range_n_N']} / {a['common_range_n_T']} | {common} |")
md+=['','Hasta 1 ve 3 negatif farkları bu duyarlılık kontrollerinde aynı yönde kalıyor. Hasta 2 ham +0.1526 iken derinlik ve döngü kontrolünde −0.0079, kimlik ve zonasyonla +0.0207 oluyor. Karışık kümeler dahil edilince +0.0727. Bu değişimler, modelin neyi çıkardığına ve hücre seçimlerine bağlılığı gösterir; kontrol değişkenlerinin nedensel açıklama olduğunu göstermez.','','Ortak aralık: beş kovaryatın iki dokudaki %5–%95 aralık kesişimi; en az 20 çekirdek/doku şartı. Yalnızca ayrı ayrı kovaryat aralıkları incelenir, çok boyutlu eşleştirme veya denge garantisi değildir. Hasta 2’de ana seçimden 3 komşu doku ve 36 tümör çekirdeği kaldığından kısıtlı model kurulmadı.','','## Hücre tipi seçimi','','Yazar kodundaki hepatosit küme numaraları, yeniden üretilmemiş özgün kümelendirmeye ait olduğundan doğrudan kullanılmadı. 2.000 gen üzerinden yeni kümelendirme yapıldı; 47 imza geni, kaynak kontrol listeleri ve MT-/RPL/RPS genleri küme oluşturma özelliklerinden çıkarıldı. Küme seçimi onarım farkları hesaplanmadan önce çoklu belirteç profilleri üzerinden kaydedildi.','','Ana kümeler: 0, 1, 3, 4, 5, 9, 10, 11, 12. Karışık 15 ve 16 yalnızca duyarlılık analizinde eklendi. 6 ve 17 belirsiz/nöral benzeri profiller nedeniyle ana analizde yok. Bunların kesin nöral hücre veya çift hücre olduğu iddia edilmez. ALB pek çok grupta bulunduğu için tek belirteç olarak kullanılmadı. PTPRC, PECAM1/VWF ve KRT19 gibi alternatif soy belirteçleri kontrol edildi. Ambient RNA ve çift hücre sorunu tamamen çözülmüş değildir.','','Katı alt küme: ana kümelerde HNF4A/ASGR1/FABP1/ALDOB genlerinden en az ikisinin saptanması ve PTPRC/PECAM1/VWF/KRT19/COL1A1 saptanmaması. Bu seçim derinlik ve dropout etkisine açıktır; yeni bir hedef hücre popülasyonu oluşturur, daha doğru anotasyon garantisi değildir.','','## Hesaplama ve kalite','','39.995 barkod metadata ile tam eşleşiyor. Matris: 58.100 özellik × 39.995 çekirdek; 54.475.478 pozitif tamsayı kayıt. Veriler önceden DIEM filtresinden geçmiş. Ek en az 200 gen ve en fazla %20 mitokondri sayımı koşulu hiçbir çekirdeği çıkarmadı. Bu, tüm teknik kalite sorunlarının çözüldüğü anlamına gelmez.','','Her örnekte log1p(sayım/toplam×10.000), gen listelerinde eşit ağırlıklı ortalama kullanıldı. CYTOR/LINC00152 ve PRAG1/SGK223 daha önce doğrulanan eş ad eşleştirmeleriyle 47 genin tamamı kapsandı. S 42, G2/M 54; kimlik 6 ve zonasyon 25 gen. Eksik MLF1IP iki kohortta da puana dahil edilmedi. Seurat AddModuleScore veya CellCycleScoring uygulanmadı.','','Her hasta için iki dokuya eşit toplam ağırlık veren doğrusal modeller kuruldu. Ana model doku göstergesi, log1p toplam sayım, S, G2/M, kimlik ve portal eksi merkezi zonasyon puanından oluşuyor. Tasarım matrisleri tam sırada; ana modellerin koşul sayıları 3.68, 7.93, 3.71. Hücre düzeyindeki katsayılar betimleyicidir; hücreler bağımsız hasta gibi test edilmedi.','','Kümelendirme için altı örneğin normalize dispersiyon sıralarıyla 2.000 özellik, merkezlenmiş/standartlaştırılmış ve ±10 kırpılmış PCA, 30 bileşen, 20-komşu grafiği, Leiden çözünürlük 1.0 ve tohum 20261010 kullanıldı. Batch düzeltmesi yapılmadı; bazı kümeler örnek ağırlıklı. Yeni anotasyon özgün yazar etiketlerinin tam tekrarı değildir.','','## Makale için çıkarım sınırı','','İki kohortta ortak, kontrollerden bağımsız artış iddiası desteklenmedi. Bulgular, hasta heterojenliği ile hücre tipi seçimi ve kontrol modeline duyarlılık açısından sunulabilir. İlk kohort scRNA-seq, ikinci kohort snRNA-seq olduğundan mutlak puanlar birleştirilmedi. İkinci kohort yalnızca üç hastadır ve yağlı karaciğer hastalığı ilişkili HCC’yi temsil eder. Hepatosit benzeri soy anotasyonu malignite kanıtı değildir. CNV/DNA ile malignite ayrımı ve daha geniş örneklem hâlâ gerekli.','','## Şekiller','','![Soy belirteçleri](GSE189175_marker_audit.png)','','![Kontrol duyarlılığı](GSE189175_control_sensitivity.png)','','[GSE189175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189175) · [Yazar analiz kodu](https://github.com/marcalva/hcc_sc_2022) · [Özgün makale](https://doi.org/10.1186/s13073-022-01055-5) · [Seurat döngü gen kaynağı](https://satijalab.org/seurat/articles/cell_cycle_vignette.html)']
Path('outputs/GSE189175_betimleyici_sinama.md').write_text('\n'.join(md)+'\n')
draft="""# Manuscript results draft — external exploratory cohort

This section reports exploratory reannotation and descriptive sensitivity analyses. It does not establish malignant-cell specificity or a causal effect.

## External cohort

We examined GSE189175, comprising 39,995 nuclei from paired tumor and adjacent non-tumor tissues from three patients with fatty-liver-disease-associated HCC. Because the deposited metadata did not provide cell-type assignments, we reannotated the dataset using unsupervised clustering and multi-lineage marker profiles. Repair signature genes and the specified control panels were excluded from clustering features. Nine hepatic-marker-supported clusters containing 25,977 nuclei were selected before calculating tissue contrasts. These assignments represent inferred lineage labels rather than reproduced author annotations or validated malignancy labels.

The mean 47-gene repair score was lower in tumor tissue in patients 1 and 3 and higher in patient 2. Unadjusted tumor-minus-adjacent differences were −0.0376, +0.1526 and −0.1472, respectively. Within-patient, equally tissue-weighted descriptive models incorporating transcript depth, S-phase and G2/M expression, hepatic identity and a portal-minus-central zonation contrast yielded corresponding tissue coefficients of −0.0271, +0.0207 and −0.0603.

Negative differences in patients 1 and 3 retained their direction under stricter marker-based selection and marginal common-range restriction. Patient 2 was sensitive to model specification: its tissue coefficient was −0.0079 with depth and cell-cycle controls and became +0.0207 after identity and zonation controls were added. Including mixed-lineage clusters increased the latter coefficient to +0.0727. Only three adjacent-tissue nuclei remained within the joint intersection of the five marginal common ranges for patient 2, preventing a restricted model under the prespecified minimum of 20 nuclei per tissue.

These observations did not support a uniform tumor-associated increase across the patients examined. They highlight sensitivity to patient context, lineage selection and covariate specification. The small independent sample size, reannotation uncertainty, unresolved ambient RNA and doublet effects, and differences between single-cell and single-nucleus assays preclude a definitive replication or mechanistic interpretation. No cell-level significance tests or causal estimates were reported.
"""
Path('outputs/manuscript_external_results_draft.md').write_text(draft)
print('Reports and figures saved.')
