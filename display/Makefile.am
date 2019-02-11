all: mad.jar

JAVA_SOURCES = \
	mad/AlignDAG.java \
	mad/Alignment.java \
	mad/AlignmentPanel.java \
	mad/Alignments.java \
	mad/BatchDocument.java \
	mad/JpegImagesToMovie.java \
	mad/MAD.java \
	mad/MadApplet.java \
	mad/MadPanel.java \
	mad/Node.java \
	mad/ProbabilityMatrices.java \
	mad/PropertyChangeHandler.java \
	mad/PropertyChangeIDs.java \
	mad/SaveAsFastaAction.java \
	mad/SaveAsTiffAction.java \
	mad/SaveAsMovAction.java \
	mad/SparseMatrix.java

# Root directory for class files (to be passed as -d option to javac)
JAVAROOT = .

# JMF directory
JMF_DIR = JMF-2.1.1e

# JAE directory
JAI_DIR = jai-1_1_3

# Options to pass to javac
AM_JAVACFLAGS = -extdirs $(srcdir)/$(JMF_DIR)/lib:$(srcdir)/$(JAI_DIR)/lib \
		-source 1.5

# Other files that should go into distribution package
EXTRA_DIST = \
	mad \
	classnoinst.stamp \
	mad.jar \
	mad/legend.jpg \
	$(JMF_DIR)/lib/customizer.jar \
	$(JMF_DIR)/lib/jmf.jar \
	$(JMF_DIR)/lib/jmf.properties \
	$(JMF_DIR)/lib/mediaplayer.jar \
	$(JMF_DIR)/lib/multiplayer.jar \
	$(JMF_DIR)/readme.html \
	$(JAI_DIR)/COPYRIGHT-jai.txt \
	$(JAI_DIR)/LICENSE-jai.txt \
	$(JAI_DIR)/DISTRIBUTIONREADME-jai.txt \
	$(JAI_DIR)/THIRDPARTYLICENSEREADME-jai.txt \
	$(JAI_DIR)/lib/jai_core.jar \
	$(JAI_DIR)/lib/jai_codec.jar

# Only specify rules for building MAD if javac is available
if HAVE_JAVAC
dist_noinst_JAVA = $(JAVA_SOURCES)

mad.jar: mad/manifest.mf \
	$(JMF_DIR)/lib/jmf.jar \
	$(JAI_DIR)/lib/jai_core.jar \
	$(JAI_DIR)/lib/jai_codec.jar \
	classnoinst.stamp
	jar cmf $< $@ \
		-C $(srcdir) mad \
		-C $(srcdir) $(JMF_DIR)/lib/jmf.jar \
		-C $(srcdir) $(JAI_DIR)/lib/jai_core.jar \
		-C $(srcdir) $(JAI_DIR)/lib/jai_codec.jar
	jar uf $@ mad

# Jar files should be cleaned as well
CLEANFILES = mad.jar
else
EXTRA_DIST += $(JAVA_SOURCES)
endif #HAVE_JAVAC

# Custom clean hook to remove class files
clean-local:
	-find . -name "*.class" -exec rm {} \;
