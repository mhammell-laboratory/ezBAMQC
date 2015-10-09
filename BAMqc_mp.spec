# -*- mode: python -*-
a = Analysis(['BAMqc_mp', 'libBAMQC/Coverage_prof.py', 'libBAMQC/GeneFeatures.py', 'libBAMQC/InerDist_prof.py', 'libBAMQC/IntervalTree.py', 'libBAMQC/Mappability.py', 'libBAMQC/parseBAM.py', 'libBAMQC/ReadDup_prof.py', 'libBAMQC/Results.py', 'libBAMQC/rRNA.py'],
             pathex=['/home/dmolik/BAM_QC/BAMQC-0.5'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='BAMqc_mp',
          debug=False,
          strip=None,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='BAMqc_mp')
