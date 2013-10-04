import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("HLTTest.MyCandidates.mycandidates_cfi")
process.isoDump.RootFileName = cms.string("hltD.root")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            #fileNames = cms.untracked.vstring('/store/caf/user/sani/ZElectrons_RunD/outputA_6_1_w2s.root'
                            fileNames = cms.untracked.vstring(
'/store/caf/user/sani/ZElectrons_RunD/outputA_10_1_a08.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_11_1_7ZI.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_12_1_TCy.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_13_1_7m3.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_14_1_Yl4.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_15_1_Czk.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_16_1_Qia.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_17_1_Mcn.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_18_1_nfr.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_19_1_hug.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_1_1_zXy.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_20_1_nwc.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_21_1_FzF.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_22_1_ofs.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_23_1_dCw.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_24_1_3c7.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_25_1_Ovq.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_26_1_47y.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_27_1_k8s.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_28_1_QT1.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_29_1_bDE.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_2_1_ilr.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_30_1_2Z1.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_31_1_Etw.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_32_1_cRP.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_33_1_Qo6.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_34_1_en5.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_35_1_Ahc.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_36_1_t1x.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_37_1_4fZ.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_38_1_bf5.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_39_1_94V.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_3_1_4Lr.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_40_1_FB7.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_41_1_Qyu.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_42_1_BlI.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_43_1_Tcz.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_44_1_JEa.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_45_1_Geu.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_46_1_hpM.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_47_1_scG.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_48_1_i0h.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_49_1_gRO.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_4_1_toZ.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_50_1_z4q.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_51_1_9to.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_52_1_CdP.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_53_1_0MI.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_54_1_nZH.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_55_1_n2G.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_56_1_Xfx.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_57_1_Hmr.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_58_1_Ku9.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_59_1_qJk.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_5_1_yuf.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_60_1_daM.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_6_1_HKD.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_7_1_M0H.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_8_1_dNK.root',
'/store/caf/user/sani/ZElectrons_RunD/outputA_9_1_kuo.root',
                                                              )
                            )

process.p = cms.Path(process.isoDump)
