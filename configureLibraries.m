function configureLibraries
% CONFIGURELIBRARIES Configure Orekit data provider
import java.io.File;
import org.orekit.data.*;

DM = DataProvidersManager.getInstance();
% TODO For linux, /data/ to be used
crawler = DirectoryCrawler(File([cd ,'\data\']));
DM.clearProviders();
DM.addProvider(crawler);

end

