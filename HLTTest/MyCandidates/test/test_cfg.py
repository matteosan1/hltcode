import FWCore.ParameterSet.Config as cms

process = cms.Process("ExREG")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V10::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_10_1_9BX.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_11_1_Qyl.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_12_1_eFI.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_13_1_0ZS.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_14_1_1kd.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_15_1_4C2.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_16_1_98k.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_17_1_BqW.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_18_1_x5V.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_19_1_WJ8.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_1_1_KcZ.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_20_1_qa7.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_21_1_4Bj.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_22_1_ndx.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_23_1_qYw.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_24_1_4fy.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_25_1_h0q.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_26_1_Wdi.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_27_1_XW3.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_28_1_u8N.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_29_1_zer.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_2_1_Sth.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_30_1_gop.root',
###'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_31_1_9Vy.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_32_1_9zP.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_33_1_ilc.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_34_1_VbC.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_35_1_hDg.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_36_1_WAE.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_37_1_wxE.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_38_1_JTZ.root',
##'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_39_1_3Z5.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_3_1_g8A.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_40_1_uvf.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_41_1_6Ku.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_42_1_akg.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_43_1_b2H.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_44_1_8da.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_45_1_CND.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_46_1_10i.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_47_1_4di.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_48_1_Qkc.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_49_1_SOw.root',
##'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_4_1_U1M.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_50_1_zKS.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_51_1_U5T.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_52_1_FPj.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_53_1_1DV.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_54_1_VbK.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_55_1_BoM.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_56_1_l2X.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_57_1_ItL.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_58_1_vrZ.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_59_1_Mww.root',
##'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_5_1_y5H.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_60_1_wKv.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_61_1_vdn.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_62_1_pl0.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_63_1_7hx.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_64_1_7eQ.root',
##'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_65_1_RJx.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_66_1_yB5.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_67_1_7zm.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_68_1_SdE.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_69_1_Iqm.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_6_1_IPE.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_70_1_xnN.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_7_1_a3a.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_8_1_y17.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_std_v4/oot_std_9_1_4sN.root',


#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_10_1_Vp6.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_11_1_mjC.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_12_1_hGs.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_13_1_HbU.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_14_1_dIT.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_15_1_G22.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_16_1_EP5.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_17_1_hUD.root',
##'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_18_1_icK.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_19_1_JqQ.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_1_1_Isy.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_20_1_JQ1.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_21_1_ceH.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_22_1_0x6.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_23_1_4ej.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_24_1_dla.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_25_1_9QS.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_26_1_dzE.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_27_1_MQQ.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_28_1_4VL.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_29_1_gNU.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_2_1_Um5.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_30_1_qjr.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_31_1_9hK.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_32_1_XEi.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_33_1_5YE.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_34_1_V7F.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_35_1_Rtf.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_36_1_kUS.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_37_1_BkL.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_38_1_AAz.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_39_1_b4U.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_3_1_0NG.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_40_1_26U.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_41_1_cp8.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_42_1_Zyv.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_43_1_tMo.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_44_1_UyG.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_45_1_qGm.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_46_1_SxY.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_47_1_q9B.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_48_1_tRD.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_49_1_Iq3.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_4_1_Nz9.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_50_1_qfi.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_51_1_Qkz.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_52_1_RHW.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_53_1_z30.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_54_1_g3C.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_55_1_L78.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_56_1_GYY.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_57_1_PCW.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_58_1_VMP.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_59_1_TxN.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_5_1_Xiy.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_60_1_P6X.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_61_1_BXP.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_62_1_IVn.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_63_1_h8D.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_64_1_TTs.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_65_1_e2X.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_66_1_OH5.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_67_1_bfm.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_68_1_Hd4.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_69_1_yAu.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_6_1_3AG.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_70_1_4Fv.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_7_1_2bE.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_8_1_KDC.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_cut_v4/oot_cut_9_1_lJY.root',

#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_10_1_STQ.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_11_1_efQ.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_12_1_Wio.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_13_1_7xn.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_15_1_ts2.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_16_1_izj.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_17_1_epR.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_18_1_CuL.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_19_1_dlm.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_1_1_SGs.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_20_1_0JL.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_21_1_Xev.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_22_1_LVc.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_23_1_2Iu.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_24_1_59S.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_25_1_AOB.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_26_1_ZKl.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_27_1_uAL.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_28_1_j8v.root',
##'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_29_1_9np.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_2_1_Jjc.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_30_1_xXa.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_31_1_eaL.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_32_1_Qld.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_33_1_g2o.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_34_1_wO5.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_35_1_z1e.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_36_1_tl7.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_37_1_mVs.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_38_1_18o.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_39_1_FvR.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_3_1_Uq2.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_40_1_0Ms.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_41_1_963.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_42_1_gqy.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_43_1_3Aj.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_44_1_I8N.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_45_1_wjF.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_46_1_trU.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_47_1_4oC.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_48_1_hDO.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_49_1_L29.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_4_1_o8o.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_50_1_k0p.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_51_1_knW.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_52_1_zhq.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_53_1_twf.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_54_1_rJx.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_55_1_hg7.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_56_1_OeO.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_57_1_5Bo.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_58_1_uWX.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_59_1_Szo.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_5_1_YXz.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_60_1_D6f.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_61_1_8NK.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_62_1_olg.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_63_1_VGn.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_64_1_Hpi.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_65_1_97p.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_66_1_eaf.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_67_1_Q9R.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_68_1_XQ5.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_69_1_gW3.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_6_1_ddV.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_70_1_YmF.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_7_1_0QI.root',
#'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_8_1_L07.root',
##'file:/hadoop/cms/store/user/sani/hlt/2014/zee_oot_multifit_v4/oot_multifit_9_1_Prw.root',

'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_10_1_Ghw.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_11_1_k2H.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_12_1_ZMp.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_13_1_eKl.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_14_1_Nq8.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_15_1_RfK.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_16_1_UXN.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_17_1_VHB.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_18_1_cHQ.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_19_1_ese.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_1_1_mVj.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_20_1_7cS.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_21_1_YSG.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_22_1_KOx.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_23_1_m1B.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_25_1_TFj.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_26_1_yFu.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_27_1_YO6.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_28_1_mGZ.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_29_1_bFj.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_2_1_jsA.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_30_1_lLi.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_31_1_wCS.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_32_1_UQv.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_33_1_yPV.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_34_1_Nav.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_35_1_WVa.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_36_1_e8l.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_37_1_uKQ.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_38_1_NnP.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_39_1_UiC.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_3_1_jjn.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_40_1_4No.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_41_1_MWb.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_42_1_b6B.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_43_1_LgG.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_44_1_eBs.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_45_1_m5e.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_46_1_3Hs.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_47_1_OKK.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_48_1_OYN.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_49_1_UHj.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_4_1_juf.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_50_1_UHw.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_51_1_Be0.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_52_1_d2K.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_53_1_63k.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_54_1_3u3.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_55_1_oOw.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_56_1_eWr.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_57_1_TZQ.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_58_1_Xku.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_59_1_a8o.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_5_1_Rbp.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_60_1_0Fh.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_61_1_Pum.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_62_1_k0h.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_63_1_z5h.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_64_1_Zfe.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_65_1_5nG.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_66_1_Jbc.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_67_1_otx.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_68_1_fKv.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_69_1_pWi.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_6_1_5MU.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_70_1_Dz3.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_7_1_TRW.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_8_1_fTo.root',
'file:/hadoop/cms/store/user/sani/hlt/2014/qcd_oot_multifit/oot_multifit_9_1_1Rc.root',

)
                            )

process.plotDistr = cms.EDAnalyzer("plotDistr",
                                   OutputFileName = cms.string("oot_qcd_multifit.root"),
                                   hitEBLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
                                   hitEELabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
                                   isData = cms.bool(False),
                                   activateNewClustering = cms.bool(True),
                                   activateOldClustering = cms.bool(False),
                                   saveReco = cms.bool(False),
                                   saveUnseeded = cms.bool(True),
                                   trgResults = cms.InputTag("TriggerResults","","TEST"),
                                   trgSelection = cms.vstring("HLT_Ele27WP85_Gsf_v1", "HLT_Ele20WP60_Ele8_Mass55_v1","HLT_Ele25WP60_SC4_Mass5_v1",)
                                   )
process.p = cms.Path(process.plotDistr)

                                        
