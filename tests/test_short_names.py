from mcnaughton2000_to_nwb.short_names import assign_short_names, parse_long_name

# Long names in the order of ANALYZED_(original_files)/FD4RAT1/TT1/WS_FTP.LOG,
# and the short names under which those files are archived.
RMAP_LONG = (
    "ES.cell1.rmap ES.cell10.rmap ES.cell11.rmap ES.cell2.rmap ES.cell3.rmap ES.cell5.rmap "
    "ES.cell6.rmap ES.cell7.rmap ES.cell8.rmap ES.cell9.rmap MC.cell1.rmap MC.cell10.rmap "
    "MC.cell11.rmap MC.cell2.rmap MC.cell4.rmap MC.cell5.rmap MC.cell6.rmap MC.cell7.rmap "
    "MC.cell8.rmap MC.cell9.rmap"
).split()
RMAP_SHORT = (
    "ES2BC0~1.RMA ES2BC4~1.RMA ES2BC6~1.RMA ES3BC8~1.RMA ES3BCA~1.RMA ES3BCC~1.RMA ESCELL~1.RMA "
    "ESCELL~2.RMA ESCELL~3.RMA ESCELL~4.RMA MC2DC0~1.RMA MC2DC2~1.RMA MC2DC4~1.RMA MC2DC6~1.RMA "
    "MC2DCE~1.RMA MC3DCC~1.RMA MCCELL~1.RMA MCCELL~2.RMA MCCELL~3.RMA MCCELL~4.RMA"
).split()

# The same for the ES and MC spike files of 1_RAW_(original_files)/FD4RAT1/TT1.
CEL_LONG = (
    "ES.cell1 ES.cell10 ES.cell11 ES.cell2 ES.cell3 ES.cell4 ES.cell5 ES.cell6 ES.cell7 ES.cell8 "
    "ES.cell9 MC.cell1 MC.cell10 MC.cell11 MC.cell2 MC.cell3 MC.cell4 MC.cell5 MC.cell6 MC.cell7 "
    "MC.cell8 MC.cell9"
).split()
CEL_SHORT = (
    "ES033D~1.CEL ES133D~1.CEL ES233D~1.CEL ES5999~1.CEL ES59A9~1.CEL ESA23D~1.CEL ESB23D~1.CEL "
    "ESC23D~1.CEL ESD23D~1.CEL ESE23D~1.CEL ESF23D~1.CEL MC034D~1.CEL MC41A9~1.CEL MC41B9~1.CEL "
    "MC824D~1.CEL MC924D~1.CEL MCA24D~1.CEL MCB24D~1.CEL MCC24D~1.CEL MCD24D~1.CEL MCE24D~1.CEL "
    "MCF24D~1.CEL"
).split()


def test_rate_map_short_names_match_archive():
    short_to_long = assign_short_names(RMAP_LONG)
    assert sorted(short_to_long) == sorted(RMAP_SHORT)
    # The first four files with a common prefix get ~1 to ~4 in order of creation.
    assert short_to_long["ESCELL~1.RMA"] == "ES.cell1.rmap"
    assert short_to_long["ESCELL~4.RMA"] == "ES.cell2.rmap"


def test_spike_file_short_names_match_archive():
    assert sorted(assign_short_names(CEL_LONG)) == sorted(CEL_SHORT)


def test_parse_long_name():
    assert parse_long_name("BL2.cell5") == ("BL2", 5)
    assert parse_long_name("ES.cell11.rmap") == ("ES", 11)
